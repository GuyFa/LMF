import math

import numpy as np
import igl
import pandas as pd

from compute_stats import compute_stats_from_V
from config import get_config, get_cols_names_values_list

import os

from utils.net_test_info_utils import TestStep3DMetricTensorOutput
from utils.run_CM_utils import (
    run_CM,
    run_CM_with_customized_init,
    run_CM_with_target_energy,
    run_CM_with_target_runtime,
)
from utils.run_GIF_utils import (
    run_GIF,
)


def check_nan_loss(target_V, source_T, pred_V):
    target_N = igl.per_vertex_normals(target_V, source_T)
    pred_N = igl.per_vertex_normals(pred_V, source_T)
    dot_N = np.sum(pred_N * target_N, 1)
    dot_N = np.clip(dot_N, 0, 1)  # adding this to avoid Nans
    angle_N = np.arccos(dot_N)
    angle_sum = np.sum(angle_N)
    array_has_nan = np.isnan(angle_sum)
    if array_has_nan:
        print("loss is nan during angle validation!")


def write_obj_with_UVs(filename, vertices, faces, uvs):
    # This function saves an obj file with UVs, with double precision for the UVs
    with open(filename, "w") as file:
        # Write vertices
        for vertex in vertices:
            file.write(f"v {' '.join(f'{coord:.16f}' for coord in vertex)}\n")

        # Write texture coordinates
        for uv_coord in uvs:
            file.write(f"vt {' '.join(f'{coord:.16f}' for coord in uv_coord)}\n")

        # Write faces
        for face in faces:
            # OBJ format uses 1-based indexing for vertices and texture coordinates
            file.write(f"f {' '.join(f'{face[i] + 1}/{face[i] + 1}' for i in range(len(face)))}\n")


def save_OBJs(sfile, tpath, source_V, source_T, target_V, pred_V):
    igl.write_obj(sfile + ".obj", source_V, source_T)

    igl.write_obj(tpath + "_target.obj", target_V, source_T)
    igl.write_obj(tpath + "_pred.obj", pred_V, source_T)
    write_obj_with_UVs(
        tpath + "_source_textured_by_uv_target.obj",
        source_V,
        source_T,
        np.delete(target_V, 2, axis=1),
    )
    write_obj_with_UVs(
        tpath + "_source_textured_by_uv_pred.obj",
        source_V,
        source_T,
        np.delete(pred_V, 2, axis=1),
    )


def compute_stats_to_summary(
    results_df,
    tpath,
    target_V,
    pred_V,
    root_dir_test,
    source_mesh_ind,
    total_loss,
    edge_compatibility_data,
    timer,
    source_T,
    source_V,
    compare_to_GIF,
    compare_with_original_CM,
    compare_with_CM_with_our_initialization,
    compare_with_CM_until_ours_distortion,
    compare_with_CM_until_ours_runtime,
    stats_only,
    pred_flips_before_fix,
    pred_flips_after_fix,
    IDT_stats,
    root_dir_cache_CM,
    root_path_CM_exe,
    root_dir_cache_GIF,
    root_path_GIF_exe,
):
    ind = results_df.shape[0]

    config = get_config()
    source_path = root_dir_test + "\\" + source_mesh_ind + ".obj"
    cur_stats = {}
    cur_stats[config.stats_net_results_df_cols.MODEL_PATH] = source_path
    cur_stats[config.stats_net_results_df_cols.OURS_TOTAL_LOSS] = float(total_loss.cpu().numpy())
    cur_stats[config.stats_net_results_df_cols.OURS_COMPATIBILITY_ABS_ERROR_MEAN] = edge_compatibility_data[0]
    cur_stats[config.stats_net_results_df_cols.OURS_COMPATIBILITY_ABS_ERROR_STD] = edge_compatibility_data[1]
    cur_stats[config.stats_net_results_df_cols.OURS_COMPATIBILITY_REL_ERROR_MEAN] = edge_compatibility_data[2]
    cur_stats[config.stats_net_results_df_cols.OURS_COMPATIBILITY_REL_ERROR_STD] = edge_compatibility_data[3]

    cur_stats[config.stats_net_results_df_cols.OURS_NETWORK_TIME_CONSUMPTION] = timer.total_time - timer.run_GIF

    cur_stats[config.stats_net_results_df_cols.OURS_GIF_TIME_CONSUMPTION] = timer.run_GIF

    cur_stats[config.stats_net_results_df_cols.OURS_TOTAL_TIME_CONSUMPTION] = timer.total_time + IDT_stats[8]

    cur_stats[config.stats_net_results_df_cols.VERTICES_COUNT] = source_V.shape[0]
    cur_stats[config.stats_net_results_df_cols.FACES_COUNT] = source_T.shape[0]
    cur_stats[config.stats_net_results_df_cols.BOUNDARY_VERTICES_COUNT] = igl.boundary_loop(source_T).shape[0]

    cur_stats[config.stats_net_results_df_cols.OURS_IDT_EDGE_FLIPS_COUNT] = IDT_stats[0]

    cur_stats[config.stats_net_results_df_cols.OURS_IDT_TIME_CONSUMPTION] = IDT_stats[8]

    (
        cur_stats[config.stats_net_results_df_cols.GT_E_SD],
        cur_stats[config.stats_net_results_df_cols.GT_K],
        cur_stats[config.stats_net_results_df_cols.GT_FLIPS],
    ) = compute_stats_from_V(source_V, source_T, target_V)

    if IDT_stats[0] == 0:
        (
            cur_stats[config.stats_net_results_df_cols.OURS_E_SD],
            cur_stats[config.stats_net_results_df_cols.OURS_K],
            _,
        ) = compute_stats_from_V(source_V, source_T, pred_V)

    else:
        cur_stats[config.stats_net_results_df_cols.OURS_E_SD] = 0
        cur_stats[config.stats_net_results_df_cols.OURS_K] = 0

    cur_stats[config.stats_net_results_df_cols.OURS_AFTER_MVC_FIX_FLIPS] = pred_flips_after_fix

    cur_stats[config.stats_net_results_df_cols.OURS_MEAN_VALUE_FIX_WAS_APPLIED] = (
        pred_flips_before_fix > 0 and IDT_stats[0] > 0
    )

    if compare_to_GIF:
        res_GIF = run_GIF(source_V, source_T, root_path_GIF_exe, root_dir_cache_GIF, interior_faces=500)
        uvs_from_GIF, pred_flips_before_fix_GIF, time_consmption_GIF = (
            res_GIF[0],
            res_GIF[1],
            res_GIF[4],
        )

        cur_stats[config.stats_net_results_df_cols.GIF_MEAN_VALUE_FIX_WAS_APPLIED] = pred_flips_before_fix_GIF > 0
        cur_stats[config.stats_net_results_df_cols.GIF_TIME_CONSUMPTION] = time_consmption_GIF
        (
            cur_stats[config.stats_net_results_df_cols.GIF_E_SD],
            cur_stats[config.stats_net_results_df_cols.GIF_K],
            cur_stats[config.stats_net_results_df_cols.GIF_AFTER_MVC_FIX_FLIPS],
        ) = compute_stats_from_V(source_V, source_T, uvs_from_GIF)

        if not stats_only:
            write_obj_with_UVs(
                tpath + "_source_textured_by_uv_GIF.obj",
                source_V,
                source_T,
                uvs_from_GIF,
            )
    else:
        cur_stats[config.stats_net_results_df_cols.GIF_E_SD] = 0
        cur_stats[config.stats_net_results_df_cols.GIF_K] = 0
        cur_stats[config.stats_net_results_df_cols.GIF_AFTER_MVC_FIX_FLIPS] = 0
        cur_stats[config.stats_net_results_df_cols.GIF_MEAN_VALUE_FIX_WAS_APPLIED] = 0
        cur_stats[config.stats_net_results_df_cols.GIF_TIME_CONSUMPTION] = 0

    if compare_with_original_CM:
        uvs_from_CM, CM_time, CM_iters, _ = run_CM(source_V, source_T, root_path_CM_exe, root_dir_cache_CM)
        cur_stats[config.stats_net_results_df_cols.CM_TIME] = CM_time
        cur_stats[config.stats_net_results_df_cols.CM_ITERS] = CM_iters
        (
            cur_stats[config.stats_net_results_df_cols.CM_E_SD],
            cur_stats[config.stats_net_results_df_cols.CM_K],
            cur_stats[config.stats_net_results_df_cols.CM_FLIPS],
        ) = compute_stats_from_V(source_V, source_T, uvs_from_CM)

        if not stats_only:
            write_obj_with_UVs(
                tpath + "_source_textured_by_uv_CM.obj",
                source_V,
                source_T,
                uvs_from_CM,
            )
    else:
        cur_stats[config.stats_net_results_df_cols.CM_TIME] = 0
        cur_stats[config.stats_net_results_df_cols.CM_ITERS] = 0
        cur_stats[config.stats_net_results_df_cols.CM_E_SD] = 0
        cur_stats[config.stats_net_results_df_cols.CM_K] = 0
        cur_stats[config.stats_net_results_df_cols.CM_FLIPS] = 0

    if compare_with_CM_with_our_initialization:
        (
            _,
            CM_time_using_init,
            CM_iters_using_init,
            e_sd_using_init,
        ) = run_CM_with_customized_init(source_V, source_T, pred_V, root_path_CM_exe, root_dir_cache_CM)
        cur_stats[config.stats_net_results_df_cols.CM_TIME_WITH_OURS_INIT] = CM_time_using_init
        cur_stats[config.stats_net_results_df_cols.CM_ITERS_WITH_OURS_INIT] = CM_iters_using_init
        cur_stats[config.stats_net_results_df_cols.CM_E_SD_WITH_OURS_INIT] = e_sd_using_init
    else:
        cur_stats[config.stats_net_results_df_cols.CM_TIME_WITH_OURS_INIT] = 0
        cur_stats[config.stats_net_results_df_cols.CM_ITERS_WITH_OURS_INIT] = 0
        cur_stats[config.stats_net_results_df_cols.CM_E_SD_WITH_OURS_INIT] = 0

    if compare_with_CM_until_ours_distortion:
        (
            _,
            CM_time_target_energy,
            CM_iters_target_energy,
            e_sd_after_reaching_target_energy,
        ) = run_CM_with_target_energy(
            source_V,
            source_T,
            cur_stats[config.stats_net_results_df_cols.OURS_E_SD],
            root_path_CM_exe,
            root_dir_cache_CM,
        )
        cur_stats[config.stats_net_results_df_cols.CM_TIME_TO_REACH_OURS_ENERGY] = CM_time_target_energy
        cur_stats[config.stats_net_results_df_cols.CM_ITERS_TO_REACH_OURS_ENERGY] = CM_iters_target_energy
        cur_stats[
            config.stats_net_results_df_cols.CM_E_SD_AFTER_REACHING_OURS_ENERGY
        ] = e_sd_after_reaching_target_energy
    else:
        cur_stats[config.stats_net_results_df_cols.CM_TIME_TO_REACH_OURS_ENERGY] = 0
        cur_stats[config.stats_net_results_df_cols.CM_ITERS_TO_REACH_OURS_ENERGY] = 0
        cur_stats[config.stats_net_results_df_cols.CM_E_SD_AFTER_REACHING_OURS_ENERGY] = 0

    if compare_with_CM_until_ours_runtime:
        (
            _,
            _,
            CM_iters_target_runtime,
            e_sd_after_reaching_target_runtime,
        ) = run_CM_with_target_runtime(
            source_V,
            source_T,
            cur_stats[config.stats_net_results_df_cols.OURS_TOTAL_TIME_CONSUMPTION],
            root_path_CM_exe,
            root_dir_cache_CM,
        )
        cur_stats[config.stats_net_results_df_cols.CM_ITERS_TO_REACH_OURS_RUNTIME] = CM_iters_target_runtime
        cur_stats[
            config.stats_net_results_df_cols.CM_E_SD_AFTER_REACHING_OURS_RUNTIME
        ] = e_sd_after_reaching_target_runtime
    else:
        cur_stats[config.stats_net_results_df_cols.CM_ITERS_TO_REACH_OURS_RUNTIME] = 0
        cur_stats[config.stats_net_results_df_cols.CM_E_SD_AFTER_REACHING_OURS_RUNTIME] = 0

    if results_df.empty:
        config = get_config()
        df_cols = get_cols_names_values_list(config.stats_net_results_df_cols)
        results_df = pd.DataFrame(columns=df_cols)

    ordered_values_for_df = [cur_stats[col] for col in results_df.columns.values.tolist()]
    results_df.loc[ind] = ordered_values_for_df
    return results_df


def test_save_results(
    batch_parts: TestStep3DMetricTensorOutput,
    results_df: pd.DataFrame,
    log_dir: str,
    root_dir_test: str,
    compare_to_GIF: bool,
    compare_with_original_CM: bool,
    compare_with_CM_with_our_initialization: bool,
    compare_with_CM_until_ours_distortion: bool,
    compare_with_CM_until_ours_runtime: bool,
    stats_only: bool,
    root_dir_cache_CM: str,
    root_path_CM_exe: str,
    root_dir_cache_GIF: str,
    root_path_GIF_exe: str,
):
    total_loss = batch_parts.total_loss
    if math.isnan(total_loss):
        print("loss is nan during validation!")
    batch_size = batch_parts.batch_size
    for i in range(batch_size):
        source_mesh_ind = batch_parts.source_ind[i]
        if batch_parts.IDT_stats[i][1] > 0:
            print(
                f"Object {source_mesh_ind} had a delta complex during IDT and could not finish the procedure, skipping... "
            )
        timer = batch_parts.timer
        sdir = os.path.join(log_dir, f"{source_mesh_ind}")
        print(f"writing source {source_mesh_ind}")

        sfile = os.path.join(sdir, f"{source_mesh_ind}")
        source_T = batch_parts.T[i]
        source_V = batch_parts.source_V[i]
        target_mesh_ind = batch_parts.target_inds[i]

        tpath = os.path.join(
            sdir,
            f"{target_mesh_ind}_from_{source_mesh_ind}",
        )
        if not os.path.exists(sdir) and not stats_only:
            try:
                os.mkdir(sdir)
            except Exception as e:
                print(f"had exception {e}, continuing to next source")
                continue

        if not stats_only:
            save_OBJs(
                sfile,
                tpath,
                source_V,
                source_T,
                batch_parts.target_V[i],
                batch_parts.pred_V[i],
            )

        results_df = compute_stats_to_summary(
            results_df,
            tpath,
            batch_parts.target_V[i],
            batch_parts.pred_V[i],
            root_dir_test,
            source_mesh_ind,
            total_loss,
            batch_parts.edge_compatibility_data[i],
            timer,
            source_T,
            source_V,
            compare_to_GIF,
            compare_with_original_CM,
            compare_with_CM_with_our_initialization,
            compare_with_CM_until_ours_distortion,
            compare_with_CM_until_ours_runtime,
            stats_only,
            batch_parts.pred_flips_before_fix[i],
            batch_parts.pred_flips_after_fix[i],
            batch_parts.IDT_stats[i],
            root_dir_cache_CM,
            root_path_CM_exe,
            root_dir_cache_GIF,
            root_path_GIF_exe,
        )

        excel_file_path = log_dir + "/results_summary.csv"
        results_df.to_csv(excel_file_path)
    return results_df.copy()
