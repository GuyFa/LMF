# defines the network and the train loop

import warnings
from typing import List

import matplotlib.pyplot as plt
import sys
import os

import DeformationEncoder
from DeformationDataset import DeformationDataset
from torch.utils.data import DataLoader
import torch
import pytorch_lightning as pl
from pytorch_lightning import seed_everything
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.callbacks import LearningRateMonitor
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from pytorch_lightning.callbacks import Callback
from pytorch_lightning.strategies import DDPStrategy

from pytorch_lightning.loggers import TensorBoardLogger

from args_from_cli import ArgsParams
from utils.net_utils import LossManager

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

USE_CUPY = True
has_gpu = "auto"
if USE_CUPY and torch.cuda.is_available():
    has_gpu = "gpu"


class MetricsCallback(Callback):
    """PyTorch Lightning metric callback."""

    def __init__(self):
        super().__init__()
        self.epochs_train_loss = LossManager(data_type="Training")
        self.epochs_validation_loss = LossManager(data_type="Validation")

    train_total_loss: List[float]
    train_matrix_loss: List[float]
    train_curvature_loss: List[float]

    def report_final_losses(self):
        self.epochs_train_loss.print_loss()
        self.epochs_validation_loss.print_loss()
        fig, axs = plt.subplots(1, 2)
        self.epochs_train_loss.plot_losses(axs[0])
        self.epochs_validation_loss.plot_losses(axs[1])
        plt.show()

    def on_validation_epoch_end(self, trainer, pl_module):
        self.epochs_validation_loss.total_loss.append(
            float(torch.stack(pl_module.steps_validation_loss.total_loss).mean())
        )
        self.epochs_validation_loss.print_loss()

    def on_train_epoch_end(self, trainer, pl_module):
        self.epochs_train_loss.total_loss.append(float(torch.stack(pl_module.steps_train_loss.total_loss).mean()))
        self.epochs_train_loss.print_loss()


def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)


def custom_collate(data):
    return data


def load_network_from_checkpoint(gen, args=None):
    from net_3d_mt import Net_3D_MT as MyNet

    model = MyNet.load_from_checkpoint(checkpoint_path=gen)
    if args is None:
        args = model.args
    # model.eval()
    # below we should add checks for any argument that cannot be changed retroactively for a loaded network
    # note we should also handle things like the encoder having different fields than what we specify
    # (e.g., loading pose params even if we don't want them) but that is too much hassle to check for something
    # that in any case breaks the loaded network

    loaded_normalization = model.get_layer_normalization_type()
    loaded_pointnet_normalization = model.get_pointnet_layer_normalization_type()
    model.args = args
    cur_normalization = model.get_layer_normalization_type()
    cur_pointnet_normalization = model.get_pointnet_layer_normalization_type()

    if cur_normalization != loaded_normalization:
        warnings.warn(
            f"args specify layer normalization of type {cur_normalization}, but network loaded from checkpoint"
            f" has {loaded_normalization} normalization"
        )
    if cur_pointnet_normalization != loaded_pointnet_normalization:
        warnings.warn(
            f"args specify pointnet layer normalization of type {cur_pointnet_normalization}, but network loaded from checkpoint"
            f" has {loaded_pointnet_normalization} normalization"
        )
    return model


def configure_loggings(log_dir_name, experiments_logs_path):
    logger = TensorBoardLogger(experiments_logs_path, name=log_dir_name)
    log_file_path = os.path.join(logger.log_dir, "logfile.txt")
    log_file_dir = os.path.dirname(log_file_path)
    if not os.path.exists(os.path.dirname(log_file_path)):
        try:
            os.mkdir(log_file_dir)
        except Exception as e:
            print(f"had exception {e}, continuing to next source")
    sys.stdout = open(log_file_path, "a")
    return logger


def create_data_loader(args, path, data_file, gen, use_dtype, dataset_type, batch_size):
    cur_dataset = DeformationDataset(
        path,
        data_file,
        gen.get_keys_to_load(True),
        gen.get_keys_to_load(False),
        use_dtype,
        dataset_type=dataset_type,
        args=args,
    )
    cur_loader = DataLoader(
        cur_dataset,
        batch_size=batch_size,
        collate_fn=custom_collate,
        pin_memory=(args.unpin_memory is None),
        drop_last=True,
        shuffle=(args.mode == "train"),
        num_workers=args.workers,
        persistent_workers=args.workers > 0,
    )
    code_length = gen.get_code_length(cur_dataset)
    point_dim = cur_dataset.get_point_dim()
    return cur_loader, code_length, point_dim


def main(gen, log_dir_name, args: ArgsParams):
    from net_3d_mt import Net_3D_MT as MyNet

    logger = configure_loggings(log_dir_name, args.experiments_logs_path)

    seed_everything(48, workers=True)

    checkpoint_callback = ModelCheckpoint(
        every_n_epochs=1,
        save_top_k=-1,  # , save_last=True
    )

    lr_monitor = LearningRateMonitor(logging_interval="step")

    termination_callback = EarlyStopping(
        monitor="val_loss",
        min_delta=0.001,
        patience=10000,
        verbose=False,
        mode="min",
    )
    cb = MetricsCallback()

    trainer = pl.Trainer(
        accelerator=has_gpu,
        devices=args.n_gpu,
        precision=args.precision,
        log_every_n_steps=200,
        max_epochs=10000,
        sync_batchnorm=args.n_gpu != 1,
        val_check_interval=1.0,
        logger=logger,
        plugins=None,
        accumulate_grad_batches=args.accumulate_grad_batches,
        num_sanity_val_steps=-1,
        enable_model_summary=False,
        enable_progress_bar=True,
        num_nodes=1,
        strategy=DDPStrategy(process_group_backend="gloo"),
        callbacks=[checkpoint_callback, lr_monitor, termination_callback, cb],
    )

    if trainer.precision == "16-true":
        use_dtype = torch.half
    elif trainer.precision == "32-true":
        use_dtype = torch.float
    elif trainer.precision == "64-true":
        use_dtype = torch.double
    else:
        raise Exception("trainer's precision is unexpected value")

    model = None
    if isinstance(gen, str):
        model = load_network_from_checkpoint(gen, args)
        gen = model.encoder
    gen.type(use_dtype)

    if args.mode == "train":
        train_loader, code_length, point_dim = create_data_loader(
            args,
            args.root_dir_train,
            args.data_file,
            gen,
            use_dtype,
            "train",
            args.batch_size_train,
        )
        valid_loader, _, _ = create_data_loader(
            args,
            args.root_dir_validation,
            args.data_file,
            gen,
            use_dtype,
            "validation",
            args.batch_size_val,
        )

    else:
        test_loader, code_length, point_dim = create_data_loader(
            args,
            args.root_dir_test,
            args.data_file,
            gen,
            use_dtype,
            "test",
            args.batch_size_test,
        )

    if model is None:
        assert isinstance(gen, DeformationEncoder.DeformationEncoder)
        model = MyNet(
            gen,
            code_length,
            point_dim=point_dim,
            args=args,
        )

    model.type(use_dtype)
    model.to(torch.device("cuda:0"))

    print(model.args)
    num_parameters = count_parameters(model)
    print("The number of parameters in the model is " + str(num_parameters))

    if args.mode == "train":
        try:
            trainer.fit(model, train_loader, valid_loader)
        except KeyboardInterrupt:
            cb.report_final_losses()
        cb.report_final_losses()
    else:
        trainer.test(model, test_loader)
        return
