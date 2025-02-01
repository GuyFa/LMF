import argparse
from typing import Optional, List

from pydantic import BaseModel


def get_arg_parser():
    """
    :return: arg parser with all relevant parametrs
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root-dir-train",
        help="Location of the root dir to get training data from",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--root-dir-validation",
        help="Location of the root dir to get validation data from",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--root-dir-test",
        help="Location of the root dir to get test data from",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--root-dir-cache-CM",
        help="Location of the root dir to for cache when running CM",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--root-path-CM-exe",
        help="Location of the exe file for CM",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--root-dir-cache-GIF",
        help="Location of the root dir to for cache when running GIF",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--root-path-GIF-exe",
        help="Location of the exe file for GIF",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--experiments-logs-path",
        help="Location of the root dir  ",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--allow-IDT-fix",
        help="If true, fix GIF result by apply GIF again after running IDT. Otherwise, apply the mean-value coordinates fix",
        action="store_true",
    )

    parser.add_argument(
        "--stats-only",
        help="Generate stats only when test, without saving any files per each model",
        action="store_true",
    )

    parser.add_argument(
        "--compare-with-original-GIF",
        help="Compare results with original GIF results",
        action="store_true",
    )

    parser.add_argument(
        "--compare-with-original-CM",
        help="Compare results with original CM",
        action="store_true",
    )

    parser.add_argument(
        "--compare-with-CM-with-our-initialization",
        help="Runs CM with ours result as initialization",
        action="store_true",
    )

    parser.add_argument(
        "--compare-with-CM-until-ours-distortion",
        help="Run CM until it reaches ours distortion",
        action="store_true",
    )

    parser.add_argument(
        "--compare-with-CM-until-ours-runtime",
        help="Run CM until it reaches ours runtime",
        action="store_true",
    )

    parser.add_argument(
        "--batch-size-train",
        help="Maximal number of sources  per batch in training, default is 1",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--batch-size-val",
        help="Maximal number of sources  per batch in validation, default is 1",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--batch-size-test",
        help="maximal number of sources  per batch in test, default is 1",
        type=int,
        default=1,
    )  # The Current code supports batch size of 1 at test

    parser.add_argument(
        "--weight-decay",
        help="Regularization weight decay default is 0",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--workers",
        help="Number of worker threads for data loading, default is 8",
        type=int,
        default=8,
    )
    parser.add_argument(
        "--data-file",
        help="If specified, use given path to load json file holding all pairs (relative to data dir" "s location)",
        type=str,
        default="data.json",
    )

    parser.add_argument("--unpin_memory", help="if specified, don" "t pin memory", action="store_true")
    parser.add_argument("--lr", help="learning rate, default is 1e-4", default=1e-4, type=float)
    parser.add_argument(
        "--precision",
        help="The precision we work in, should be 16,32,or 64 (default)",
        default=64,
        type=int,
    )

    parser.add_argument(
        "--interior-faces-in-GIF",
        help="The interior faces count used in GIF (only for ours). Insert -1 to use the default formula for interior faces count",
        default=-1,
        type=int,
    )

    parser.add_argument(
        "--boundary-segment-size-in-GIF",
        help="The boundary segment size used in GIF (only for ours). Insert -1 to use the default formula for boundary segment size",
        default=-1,
        type=int,
    )

    parser.add_argument(
        "--curvature-meta-vertices-rate-in-GIF",
        help="The curvature meta vertices rate used in GIF (only for ours)",
        default=0.001,
        type=float,
    )

    parser.add_argument(
        "--outer-termination-condition-rate-in-GIF",
        help="The outer termination condition rate used in GIF (only for ours)",
        default=0.1,
        type=float,
    )

    parser.add_argument(
        "--energy-related-termination-condition-in-GIF",
        help="The energy related termination condition used in GIF (only for ours)",
        default=2.0,
        type=float,
    )

    parser.add_argument("--n-gpu", help="num of gpus, default is all", type=int, default=-1)
    parser.add_argument(
        "--pointnet-layer-normalization",
        help="Type of normalization for the PointNet encoder's layers",
        default="GROUPNORM",
        type=str,
        choices={"GROUPNORM", "BATCHNORM", "IDENTITY"},
    )
    parser.add_argument(
        "--layer-normalization",
        help="Type of normalization for the decoder's layers",
        default="GROUPNORM",
        type=str,
        choices={"GROUPNORM", "GROUPNORM2", "BATCHNORM", "IDENTITY", "LAYERNORM"},
    )

    parser.add_argument("--mode", help="run option", choices=["train", "test"])

    parser.add_argument(
        "--lr-epoch-steps",
        nargs="+",
        type=int,
        help="List of epochs in which we need to decrease the learning rate",
        default=[30, 40],
    )
    parser.add_argument(
        "--accumulate-grad-batches",
        type=int,
        help="Configure the parameter for the network",
        default=1,
    )
    parser.add_argument(
        "--shuffle-triangles",
        type=bool,
        help="Shuffle triangles before NJF decoder, to avoid that group-norm overfits to a particular triangulation.",
        default=1,
    )

    parser.add_argument("--no-rotational-alignment", help="Disable the rotational alignment", action="store_true")

    return parser


class ArgsParams(BaseModel):
    root_dir_train: Optional[str]
    root_dir_validation: Optional[str]
    root_dir_test: Optional[str]
    root_dir_cache_CM: Optional[str]
    root_path_CM_exe: Optional[str]
    root_dir_cache_GIF: Optional[str]
    root_path_GIF_exe: Optional[str]
    experiments_logs_path: Optional[str]
    allow_IDT_fix: bool
    stats_only: bool
    compare_with_original_GIF: bool
    compare_with_original_CM: bool
    compare_with_CM_with_our_initialization: bool
    compare_with_CM_until_ours_distortion: bool
    compare_with_CM_until_ours_runtime: bool
    batch_size_train: int
    batch_size_val: int
    batch_size_test: int
    weight_decay: float
    workers: int
    data_file: str
    unpin_memory: bool
    lr: float
    precision: int
    interior_faces_in_GIF: int
    boundary_segment_size_in_GIF: int
    curvature_meta_vertices_rate_in_GIF: float
    outer_termination_condition_rate_in_GIF: float
    energy_related_termination_condition_in_GIF: float
    n_gpu: int
    pointnet_layer_normalization: str
    layer_normalization: str
    mode: str
    lr_epoch_steps: List[int]
    accumulate_grad_batches: int
    shuffle_triangles: bool
    no_rotational_alignment: bool


def parse_args():
    """
    parse command line args and return them
    :return: the args in stadnard arg parser's args format
    """
    parser = get_arg_parser()
    args = parser.parse_args()

    args_params = ArgsParams(**vars(args))
    assert args.interior_faces_in_GIF >= -1
    assert args.boundary_segment_size_in_GIF >= -1
    assert args.curvature_meta_vertices_rate_in_GIF >= 0
    assert args.outer_termination_condition_rate_in_GIF >= 0
    assert args.energy_related_termination_condition_in_GIF >= 0

    return args_params
