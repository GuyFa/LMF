import sys
import os
from source_lmf.Experiment import Experiment
import torch
import gc
import cupy


os.environ["PL_TORCH_DISTRIBUTED_BACKEND"] = "gloo"

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "source_lmf"))


class DefaultExperiment(Experiment):
    def __init__(self):
        Experiment.__init__(
            self,
            "DefaultExperiment",
            "This experiement takes a folder of meshes as input and a list of pairs (source_mesh, target_mesh) in this folder (in pairs.json) and learns to deform the source into the target.",
        )

    def init_encoder(self, encoder, args):
        # source|target=True indicates that the source will ba passed to a PointNet encoder
        encoder.add_pointnet(1000, source=True, target=False)


if __name__ == "__main__":
    gc.collect()
    torch.cuda.empty_cache()
    mempool = cupy.get_default_memory_pool()
    pinned_mempool = cupy.get_default_pinned_memory_pool()
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()

    mode = 1  # 0 means train, 1 means test
    if mode == 0:
        # train commands:
        sys.argv.append("--mode")
        sys.argv.append("train")

        sys.argv.append("--root-dir-train")
        sys.argv.append("C:/Data/Projects/DLGP/Datasets/train patches")

        sys.argv.append("--root-dir-validation")
        sys.argv.append("C:/Data/Projects/DLGP/Datasets/validation patches")

        sys.argv.append("--batch-size-train")
        sys.argv.append("1")
        sys.argv.append("--batch-size-val")
        sys.argv.append("10")

        sys.argv.append("--accumulate-grad-batches")
        sys.argv.append("1")

        sys.argv.append("--lr")
        sys.argv.append("1e-4")


    elif mode == 1:
        # test commands:
        sys.argv.append("--mode")
        sys.argv.append("test")

        sys.argv.append("--root-dir-test")
        sys.argv.append("C:/Data/Projects/DLGP/Datasets/Low GT energies set")

        sys.argv.append("--batch-size-test")
        sys.argv.append("1")
        # sys.argv.append("--compare-with-original-GIF")
        # sys.argv.append("--compare-with-original-CM")
        # sys.argv.append("--compare-with-CM-with-our-initialization")
        # sys.argv.append("--compare-with-CM-until-ours-distortion")
        # sys.argv.append("--compare-with-CM-until-ours-runtime")

        sys.argv.append("--root-dir-cache-CM")
        sys.argv.append("C:/Data/DL_experiments_data/DLGP/cache_CM")
        sys.argv.append("--root-path-CM-exe")
        sys.argv.append("C:/Data/Repositories/CompMajor/Parameterization_cmd.exe")

        sys.argv.append("--root-dir-cache-GIF")
        sys.argv.append("C:/Data/DL_experiments_data/DLGP/cache_GIF")
        sys.argv.append("--root-path-GIF-exe")
        sys.argv.append("C:/Data/Repositories/DLGP_GIF/GIF/x64/Release/GIF.exe")

        # sys.argv.append("--interior-faces-in-GIF")
        # sys.argv.append("500")
        # sys.argv.append("--boundary-segment-size-in-GIF")
        # sys.argv.append("1")

        sys.argv.append("--stats-only")
        # sys.argv.append("--allow-IDT-fix")

    # gpu settings:

    sys.argv.append("--n-gpu")
    sys.argv.append("1")

    sys.argv.append("--precision")
    sys.argv.append("32")

    sys.argv.append("--experiments-logs-path")
    sys.argv.append("C:/Data/DL_experiments_data/DLGP/lightning_logs")

    exp = DefaultExperiment()
    # this parses the command line arguments and then trains on the given dataset with the given experiment
    exp.get_args_and_train()
