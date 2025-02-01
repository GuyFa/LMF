import os
from datetime import datetime

import igl
from subprocess import Popen, PIPE, TimeoutExpired


TIMEOUT_SECONDS = 3600


def save_data_for_CM(vertices, faces, cache_path, init_UVs=None):
    cur_time = datetime.now().strftime("%Y%m%d%H%M%S")
    tmp_input_obj_file = cache_path + "/" + cur_time + ".obj"
    igl.write_obj(tmp_input_obj_file, vertices, faces)
    tmp_result_obj_file = cache_path + "/" + cur_time + "_CM_res.obj"
    tmp_meta_data_txt_file = cache_path + "/" + cur_time + "_CM_res_meta_data.txt"
    tmp_init_obj_file = None
    if init_UVs is not None:
        tmp_init_obj_file = cache_path + "/" + cur_time + "_init.obj"
        igl.write_obj(tmp_init_obj_file, init_UVs, faces)
    return (
        tmp_input_obj_file,
        tmp_result_obj_file,
        tmp_meta_data_txt_file,
        tmp_init_obj_file,
    )


def run_CM(V, F, exe_path, cache_path):
    (
        tmp_input_obj_file,
        tmp_result_obj_file,
        tmp_meta_data_txt_file,
        _,
    ) = save_data_for_CM(V, F, cache_path)

    exe_file_dir, exe_file_name = exe_path.rsplit("/", 1)

    cmdline = (
        exe_file_name + ' "' + tmp_input_obj_file + '" "' + tmp_result_obj_file + '" "' + tmp_meta_data_txt_file + '"'
    )

    p = Popen(
        "cmd /c " + cmdline,
        cwd=exe_file_dir,
        shell=True,
        stdout=PIPE,
        stderr=PIPE,
    )
    try:
        stdout, stderr = p.communicate(timeout=TIMEOUT_SECONDS)
    except TimeoutExpired:
        # Handle timeout
        print("Command timed out. Terminating the process.")
        p.terminate()  # Terminate the process

    uvs_full, _, _, _, _, _ = igl.read_obj(tmp_result_obj_file)
    meta_data = open(tmp_meta_data_txt_file, "r").readlines()
    CM_time_consumption = float(meta_data[0].removesuffix("\n").split(":")[1])
    CM_iters_count = float(meta_data[1].removesuffix("\n").split(":")[1])
    e_sd = float(meta_data[2].removesuffix("\n").split(":")[1])

    os.remove(tmp_input_obj_file)
    os.remove(tmp_result_obj_file)
    os.remove(tmp_meta_data_txt_file)
    return uvs_full[:, :2], CM_time_consumption, CM_iters_count, e_sd


def run_CM_with_target_energy(V, F, target_energy, exe_path, cache_path):
    (
        tmp_input_obj_file,
        tmp_result_obj_file,
        tmp_meta_data_txt_file,
        _,
    ) = save_data_for_CM(V, F, cache_path)

    exe_file_dir, exe_file_name = exe_path.rsplit("/", 1)

    cmdline = (
        exe_file_name
        + ' "'
        + tmp_input_obj_file
        + '" "'
        + tmp_result_obj_file
        + '" "'
        + tmp_meta_data_txt_file
        + '" '
        + str(target_energy)
    )

    p = Popen(
        "cmd /c " + cmdline,
        cwd=exe_file_dir,
        shell=True,
        stdout=PIPE,
        stderr=PIPE,
    )
    try:
        stdout, stderr = p.communicate(timeout=TIMEOUT_SECONDS)
    except TimeoutExpired:
        # Handle timeout
        print("Command timed out. Terminating the process.")
        p.terminate()  # Terminate the process

    uvs_full, _, _, _, _, _ = igl.read_obj(tmp_result_obj_file)
    meta_data = open(tmp_meta_data_txt_file, "r").readlines()
    CM_time_consumption = float(meta_data[0].removesuffix("\n").split(":")[1])
    CM_iters_count = float(meta_data[1].removesuffix("\n").split(":")[1])
    e_sd = float(meta_data[2].removesuffix("\n").split(":")[1])

    os.remove(tmp_input_obj_file)
    os.remove(tmp_result_obj_file)
    os.remove(tmp_meta_data_txt_file)
    return uvs_full[:, :2], CM_time_consumption, CM_iters_count, e_sd


def run_CM_with_target_runtime(V, F, target_time, exe_path, cache_path):
    (
        tmp_input_obj_file,
        tmp_result_obj_file,
        tmp_meta_data_txt_file,
        _,
    ) = save_data_for_CM(V, F, cache_path)

    exe_file_dir, exe_file_name = exe_path.rsplit("/", 1)

    cmdline = (
        exe_file_name
        + ' "'
        + tmp_input_obj_file
        + '" "'
        + tmp_result_obj_file
        + '" "'
        + tmp_meta_data_txt_file
        + '" 2 '
        + str(target_time)
    )

    p = Popen(
        "cmd /c " + cmdline,
        cwd=exe_file_dir,
        shell=True,
        stdout=PIPE,
        stderr=PIPE,
    )
    try:
        stdout, stderr = p.communicate(timeout=TIMEOUT_SECONDS)
    except TimeoutExpired:
        # Handle timeout
        print("Command timed out. Terminating the process.")
        p.terminate()  # Terminate the process

    uvs_full, _, _, _, _, _ = igl.read_obj(tmp_result_obj_file)
    meta_data = open(tmp_meta_data_txt_file, "r").readlines()
    CM_time_consumption = float(meta_data[0].removesuffix("\n").split(":")[1])
    CM_iters_count = float(meta_data[1].removesuffix("\n").split(":")[1])
    e_sd = float(meta_data[2].removesuffix("\n").split(":")[1])

    os.remove(tmp_input_obj_file)
    os.remove(tmp_result_obj_file)
    os.remove(tmp_meta_data_txt_file)
    return uvs_full[:, :2], CM_time_consumption, CM_iters_count, e_sd


def run_CM_with_customized_init(V, F, init_UVs, exe_path, cache_path):
    (
        tmp_input_obj_file,
        tmp_result_obj_file,
        tmp_meta_data_txt_file,
        tmp_init_obj_file,
    ) = save_data_for_CM(V, F, cache_path, init_UVs)

    exe_file_dir, exe_file_name = exe_path.rsplit("/", 1)

    cmdline = (
        exe_file_name
        + ' "'
        + tmp_input_obj_file
        + '" "'
        + tmp_result_obj_file
        + '" "'
        + tmp_meta_data_txt_file
        + '" 2 100000 "'
        + tmp_init_obj_file
        + '"'
    )

    p = Popen(
        "cmd /c " + cmdline,
        cwd=exe_file_dir,
        shell=True,
        stdout=PIPE,
        stderr=PIPE,
    )
    try:
        stdout, stderr = p.communicate(timeout=TIMEOUT_SECONDS)
    except TimeoutExpired:
        # Handle timeout
        print("Command timed out. Terminating the process.")
        p.terminate()  # Terminate the process

    uvs_full, _, _, _, _, _ = igl.read_obj(tmp_result_obj_file)
    meta_data = open(tmp_meta_data_txt_file, "r").readlines()
    CM_time_consumption = float(meta_data[0].removesuffix("\n").split(":")[1])
    CM_iters_count = float(meta_data[1].removesuffix("\n").split(":")[1])
    e_sd = float(meta_data[2].removesuffix("\n").split(":")[1])

    os.remove(tmp_input_obj_file)
    os.remove(tmp_result_obj_file)
    os.remove(tmp_meta_data_txt_file)
    os.remove(tmp_init_obj_file)
    return uvs_full[:, :2], CM_time_consumption, CM_iters_count, e_sd
