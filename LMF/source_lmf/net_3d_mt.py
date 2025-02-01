import numpy as np

from torch import nn
import torch
import PerCentroidBatchMaker
import time
import pytorch_lightning as pl

import math
import pandas as pd

from args_from_cli import ArgsParams

from utils.net_evaluation_utils import (
    weighted_mse_loss,
    compute_matrix_loss,
    log_mat,
    exp_mat,
)
from utils.net_prediction_utils import (
    get_M_list_from_raw_input,
)
from utils.net_test_utils import test_save_results
from utils.net_test_info_utils import TestStep3DMetricTensorOutput, Timinigs
from utils.net_utils import release_memory, LossManager
from utils.run_GIF_utils import compute_UVs_from_metric_tensors_GIF


class Net_3D_MT(pl.LightningModule):
    """
    the network
    """

    def __init__(self, encoder, code_dim, args: ArgsParams, point_dim=6):
        print("********** Some Network info...")
        print(f"********** code dim: {code_dim}")
        print(f"********** centroid dim: {point_dim}")
        super().__init__()
        self.args = args

        layer_normalization = self.get_layer_normalization_type()
        if layer_normalization == "IDENTITY":
            # print("Using IDENTITY (no normalization) in per_face_decoder!")
            self.per_face_decoder = nn.Sequential(
                nn.Linear(point_dim + code_dim, 128),
                nn.Identity(),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.Identity(),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.Identity(),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.Identity(),
                nn.ReLU(),
                nn.Linear(128, 6),
            )
        elif layer_normalization == "BATCHNORM":
            # print("Using BATCHNORM in per_face_decoder!")
            self.per_face_decoder = nn.Sequential(
                nn.Linear(point_dim + code_dim, 128),
                nn.BatchNorm1d(128),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.BatchNorm1d(128),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.BatchNorm1d(128),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.BatchNorm1d(128),
                nn.ReLU(),
                nn.Linear(128, 6),
            )
        elif layer_normalization == "GROUPNORM_CONV":
            # print("Using GROUPNORM2 in per_face_decoder!")
            self.per_face_decoder = nn.Sequential(
                nn.Conv1d(point_dim + code_dim, 128, 1),
                nn.GroupNorm(num_groups=4, num_channels=128),
                nn.ReLU(),
                nn.Conv1d(128, 128, 1),
                nn.GroupNorm(num_groups=4, num_channels=128),
                nn.ReLU(),
                nn.Conv1d(128, 128, 1),
                nn.GroupNorm(num_groups=4, num_channels=128),
                nn.ReLU(),
                nn.Conv1d(128, 128, 1),
                nn.GroupNorm(num_groups=4, num_channels=128),
                nn.ReLU(),
                nn.Conv1d(128, 6, 1),
            )
        elif layer_normalization == "GROUPNORM":
            print("Using GROUPNORM in per_face_decoder!")
            self.per_face_decoder = nn.Sequential(
                nn.Linear(point_dim + code_dim, 128),
                nn.GroupNorm(num_groups=4, num_channels=128),
                # , eps=0.0001 I have considered increasing this value in case we have channels from pointnet with the same values.
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.GroupNorm(num_groups=4, num_channels=128),
                # , eps=0.0001 I have considered increasing this value in case we have channels from pointnet with the same values.
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.GroupNorm(num_groups=4, num_channels=128),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.GroupNorm(num_groups=4, num_channels=128),
                nn.ReLU(),
                nn.Linear(128, 6),
            )
        elif layer_normalization == "LAYERNORM":
            # print("Using LAYERNORM in per_face_decoder!")
            self.per_face_decoder = nn.Sequential(
                nn.Linear(point_dim + code_dim, 128),
                nn.LayerNorm(normalized_shape=128),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.LayerNorm(normalized_shape=128),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.LayerNorm(normalized_shape=128),
                nn.ReLU(),
                nn.Linear(128, 128),
                nn.LayerNorm(normalized_shape=128),
                nn.ReLU(),
                nn.Linear(128, 6),
            )
        else:
            raise Exception("unknown normalization method")

        self.__IDENTITY_INIT = True
        if self.__IDENTITY_INIT:
            self.per_face_decoder._modules["12"].bias.data.zero_()
            self.per_face_decoder._modules["12"].weight.data.zero_()

        self.encoder = encoder
        self.point_dim = point_dim
        self.code_dim = code_dim
        self.weighted_mse_loss = weighted_mse_loss
        self.save_hyperparameters()
        self.log_validate = True
        self.lr = args.lr
        self.val_step_iter = 0
        self.df = pd.DataFrame()
        self.steps_train_loss = LossManager()
        self.steps_validation_loss = LossManager()

    ##################
    # inference code below
    ##################
    def forward(self, x):
        """
        The MLP applied to a (batch) of global code concatenated to a centroid (z|c)
        :param x: B x (|z|+|c|) batch of (z|c) vectors
        :return: B x 9 batch of 9 values that are the 3x3 matrix predictions for each input vector
        """
        if self.code_dim + self.point_dim < x.shape[1]:
            print("WARNING: discarding part of the latent code.")
            x = x[:, : self.code_dim + self.point_dim]

        return self.per_face_decoder(x.type(self.per_face_decoder[0].bias.type()))

    def predict_metric_tensors(self, source):
        """
        given a batch class, predict jacobians
        :param single_source_batch: batch object
        :return: BxTx3x3 a tensor of 3x3 jacobians, per T tris, per B targets in batch
        """
        # extract the encoding of the source and target
        codes = self.extract_code(source)
        # get the network predictions, a BxTx3x3 tensor of 3x3 metric tensors, per T tri, per B target in batch
        return self.predict_metric_tensors_from_codes(codes, source)

    def predict_metric_tensors_from_codes(self, codes, source):
        """
        predict jacobians w.r.t give global codes and the batch
        :param codes: codes for each source/target in batch
        :param single_source_batch: the batch
        :return:BxTx3x3 a tensor of 3x3 metric tensors, per T tris, per B targets in batch
        """
        # take all encodings z_i of targets, and all centroids c_j of triangles, and create a cartesian product of the two as a 2D tensor so each sample in it is a vector with rows (z_i|c_j)

        centroids_and_normals = [s.get_centroids_and_normals() for s in source]
        data = torch.cat(centroids_and_normals, dim=0)
        net_input = PerCentroidBatchMaker.PerCentroidBatchMaker(codes, data, args=self.args)
        # stacked = net_input.to_stacked()
        stacked = net_input.to_stacked_batch()
        if self.args.layer_normalization != "GROUPNORM2":
            stacked = net_input.prep_for_linear_layer(stacked)
        else:
            stacked = net_input.prep_for_conv1d(stacked)
        # feed the 2D tensor through the network, and get a 3x3 matrix for each (z_i|c_j)
        res = self.forward(stacked)
        # because of stacking the result is a 6-entry vec for each (z_i|c_j), now let's turn it to a batch x tris x 6 tensor
        pred_M = net_input.back_to_non_stacked(res)

        ret = pred_M.reshape(pred_M.shape[0], pred_M.shape[1], 6)
        if torch.isinf(ret).any():
            print("net's output contains inf!")
        if torch.isnan(ret).any():
            print("net's output contains nan!")
        if torch.isinf(ret).any() or torch.isnan(ret).any():
            return torch.nan_to_num(ret)
        else:
            return ret

    def extract_code(self, source):
        """
        given a batch, extract the global code w.r.t the source and targets, using the set encoders
        :param batch: the batch object
        :return: Bx|z| batch of codes z
        """
        codes = self.encoder.encode_deformation(source)
        faces_counts = [s.mesh_processor.faces.shape[0] for s in source]
        codes_per_traingles = torch.repeat_interleave(codes, torch.tensor(faces_counts).to(self.device), dim=0)
        return codes_per_traingles

    #######################################
    # Pytorch Lightning Boilerplate code (training, logging, etc.)
    #######################################

    def training_step(self, source_batches, batch_id):
        cur_loss = self.my_step(source_batches, batch_id)
        optimizer = self.optimizers()
        for param_group in optimizer.param_groups:
            print("learning rate: " + str(param_group["lr"]))
        self.steps_train_loss.total_loss.append(cur_loss)
        return cur_loss

    def validation_step(self, batch, batch_idx):
        cur_loss = self.my_step(batch, batch_idx)
        self.val_step_iter += 1
        # these next few lines make sure cupy releases all memory
        release_memory()
        if math.isnan(cur_loss):
            print("loss is nan during validation!")
        self.steps_validation_loss.total_loss.append(cur_loss)
        return cur_loss

    def test_step(self, batch, batch_idx):
        batch_parts = self.my_step(batch, batch_idx)
        self.df = test_save_results(
            batch_parts,
            self.df.copy(),
            self.logger.log_dir,
            self.args.root_dir_test,
            self.args.compare_with_original_GIF,
            self.args.compare_with_original_CM,
            self.args.compare_with_CM_with_our_initialization,
            self.args.compare_with_CM_until_ours_distortion,
            self.args.compare_with_CM_until_ours_runtime,
            self.args.stats_only,
            self.args.root_dir_cache_CM,
            self.args.root_path_CM_exe,
            self.args.root_dir_cache_GIF,
            self.args.root_path_GIF_exe,
        )

    def get_gt_map(self, source, target):
        f = source.mesh_processor.faces

        s_v = source.get_vertices()
        t_v = target.get_vertices()

        f_s_v = torch.zeros((f.shape[0], 3, 3), device="cuda:0")
        f_s_v[:, 0, :] = s_v[f[:, 0]]
        f_s_v[:, 1, :] = s_v[f[:, 1]]
        f_s_v[:, 2, :] = s_v[f[:, 2]]
        s_v_cross = torch.cross((f_s_v[:, 1, :] - f_s_v[:, 0, :]), (f_s_v[:, 2, :] - f_s_v[:, 0, :]))

        f_s_v_4 = f_s_v[:, 0, :] + s_v_cross / torch.sqrt(torch.linalg.norm(s_v_cross, axis=1)).reshape((-1, 1)).repeat(
            1, 3
        )
        s_j = torch.zeros((f.shape[0], 3, 3), device="cuda:0")
        s_j[:, :, 0] = f_s_v[:, 1, :] - f_s_v[:, 0, :]
        s_j[:, :, 1] = f_s_v[:, 2, :] - f_s_v[:, 0, :]
        s_j[:, :, 2] = f_s_v_4 - f_s_v[:, 0, :]

        t_v = t_v.squeeze()

        f_t_v = torch.zeros((f.shape[0], 3, 3), device="cuda:0")
        f_t_v[:, 0, :] = t_v[f[:, 0]]
        f_t_v[:, 1, :] = t_v[f[:, 1]]
        f_t_v[:, 2, :] = t_v[f[:, 2]]
        t_v_cross = torch.cross((f_t_v[:, 1, :] - f_t_v[:, 0, :]), (f_t_v[:, 2, :] - f_t_v[:, 0, :]))

        f_t_v_4 = f_t_v[:, 0, :] + t_v_cross / torch.sqrt(torch.linalg.norm(t_v_cross, axis=1)).reshape((-1, 1)).repeat(
            1, 3
        )
        t_j = torch.zeros((f.shape[0], 3, 3), device="cuda:0")
        t_j[:, :, 0] = f_t_v[:, 1, :] - f_t_v[:, 0, :]
        t_j[:, :, 1] = f_t_v[:, 2, :] - f_t_v[:, 0, :]
        t_j[:, :, 2] = f_t_v_4 - f_t_v[:, 0, :]

        inverse_s_j = torch.linalg.inv(s_j)
        j = torch.einsum("abc,acd->abd", t_j, inverse_s_j)
        m = torch.einsum("abc,acd->abd", torch.transpose(j, 2, 1), j)
        log_m = log_mat(m)
        return log_m

    def predict_map(self, source):
        faces_counts = [s.mesh_processor.faces.shape[0] for s in source]

        timer = Timinigs()
        start = time.time()
        pred_log_M_raw = self.predict_metric_tensors(source)

        pred_log_M_list = get_M_list_from_raw_input(pred_log_M_raw, faces_counts)

        end = time.time()
        timer.predict_mt = end - start

        if self.args.mode == "train":
            return None, pred_log_M_list, None

        pred_M_list = []

        start = time.time()
        for pred_log_M in pred_log_M_list:
            pred_M_list.append(exp_mat(pred_log_M.squeeze()).unsqueeze(0))
        end = time.time()
        timer.compute_exp_mt = end - start

        pred_V_list = [
            compute_UVs_from_metric_tensors_GIF(
                pred_M_list[i].squeeze(),
                source[i].mesh_processor.vertices,
                source[i].mesh_processor.faces,
                self.args.root_path_GIF_exe,
                self.args.root_dir_cache_GIF,
                self.args.interior_faces_in_GIF,
                self.args.boundary_segment_size_in_GIF,
                self.args.curvature_meta_vertices_rate_in_GIF,
                self.args.outer_termination_condition_rate_in_GIF,
                self.args.energy_related_termination_condition_in_GIF,
                self.args.allow_IDT_fix,
                False,
            )
            for i in range(len(source))
        ]

        timer.run_GIF = np.sum(np.array([pred_V_list[i][4] for i in range(len(source))]))
        timer.object_alignment_staff = source[0].mesh_processor.timings.object_alignment_staff
        timer.all_samples = source[0].mesh_processor.timings.all_samples
        timer.points_and_normals = source[0].mesh_processor.timings.points_and_normals
        timer.move_data_to_gpu = source[0].mesh_processor.timings.move_data_to_gpu

        return (
            timer,
            pred_log_M_list,
            pred_V_list,
        )

    def my_step(self, source_batch, batch_idx, test=False):
        source = [s[0] for s in source_batch]
        target = [s[1] for s in source_batch]

        (
            timer,
            pred_log_M_list,
            pred_V_list,
        ) = self.predict_map(source)

        log_m = [self.get_gt_map(source[i], target[i]) for i in range(len(source))]

        # compute total_loss
        total_loss = compute_matrix_loss(source, log_m, pred_log_M_list)

        if self.args.mode == "test":
            # Retrieve the original mesh before alignment
            for i in range(len(source)):
                source[i].mesh_processor.get_vertices_before_alignment()
                source[i].mesh_processor.get_alignment_scale()
                target[i].mesh_processor.get_vertices_before_alignment()

                start_time = time.time()
                pred_V_list[i][0] *= 1 / source[i].mesh_processor.alignment_scale
                end_time = time.time()
                timer.object_alignment_staff += end_time - start_time

            ret = TestStep3DMetricTensorOutput(
                target_V=[target[i].mesh_processor.original_vertices for i in range(len(source))],
                source_V=[source[i].mesh_processor.original_vertices for i in range(len(source))],
                total_loss=total_loss,
                T=[source[i].get_source_triangles().squeeze() for i in range(len(source))],
                source_ind=[source[i].source_ind for i in range(len(source))],
                target_inds=[target[i].target_ind for i in range(len(source))],
                batch_size=len(source),
                pred_V=[
                    np.hstack([pred_V_list[i][0], np.zeros((pred_V_list[i][0].shape[0], 1))])
                    for i in range(len(source))
                ],
                pred_flips_before_fix=np.array([pred_V_list[i][1] for i in range(len(source))]),
                pred_flips_after_fix=np.array([pred_V_list[i][2] for i in range(len(source))]),
                IDT_stats=[pred_V_list[i][3] for i in range(len(source))],
                timer=timer,
                edge_compatibility_data=[pred_V_list[i][5] for i in range(len(source))],
            )
        else:
            ret = total_loss
        return ret

    def on_validation_epoch_end(self):
        losses = torch.tensor(self.steps_validation_loss.total_loss)
        losses_arr = torch.tensor(losses, device="cpu").numpy()
        val_loss = np.mean(losses_arr)
        self.log(
            "val_loss",
            val_loss,
            batch_size=self.args.batch_size_val,
            logger=True,
            on_epoch=True,
        )
        self.val_step_iter = 0
        self.log_validate = True
        self.steps_validation_loss.total_loss.clear()

    def on_train_epoch_end(self) -> None:
        losses = torch.tensor(self.steps_train_loss.total_loss)
        losses_arr = torch.tensor(losses, device="cpu").numpy()
        train_loss = np.mean(losses_arr)
        self.log(
            "train_loss",
            train_loss,
            batch_size=self.args.batch_size_train,
            logger=True,
            on_epoch=True,
        )
        release_memory()
        self.steps_train_loss.total_loss.clear()

    def on_train_batch_end(self, outputs, batch, batch_idx, unused=0):
        release_memory()
        return

    def on_validation_batch_end(self, outputs, batch, batch_idx, unused=0):
        release_memory()
        return

    def on_test_batch_end(self, outputs, batch, batch_idx, unused=0):
        release_memory()
        return

    def get_layer_normalization_type(self):
        if hasattr(self.args, "layer_normalization"):
            layer_normalization = self.args.layer_normalization
        else:
            assert hasattr(self.args, "batchnorm_decoder")
            layer_normalization = self.args.batchnorm_decoder
        return layer_normalization

    def get_pointnet_layer_normalization_type(self):
        if hasattr(self.args, "pointnet_layer_normalization"):
            layer_normalization = self.args.pointnet_layer_normalization
        else:
            assert hasattr(self.args, "batchnorm_encoder")
            layer_normalization = self.args.batchnorm_encoder
        return layer_normalization

    def configure_optimizers(self):
        if self.args.weight_decay == 0:
            print("No weight decay regulaization!")
            optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        else:
            optimizer = torch.optim.Adam(self.parameters(), lr=self.lr, weight_decay=self.args.weight_decay)
        lr_scheduler = {
            "scheduler": torch.optim.lr_scheduler.MultiStepLR(
                optimizer,
                milestones=self.args.lr_epoch_steps,
                gamma=0.7,
            ),
            "name": "scheduler",
        }
        return [optimizer], [lr_scheduler]
