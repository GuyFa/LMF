import numpy as np
import torch


def get_M_list_from_raw_input(pred_M_raw, faces_counts):
    pred_M = torch.zeros(
        (pred_M_raw.shape[0], pred_M_raw.shape[1], 3, 3),
        device="cuda:0",
        dtype=torch.float64,
    )

    # Denote each row as (A, B, C, D, E, F). The metric tensor is
    # [[A, B, C],
    # [B, D, E],
    # [C, E, F]]
    pred_M[:, :, 0, 0] = pred_M_raw[:, :, 0]  # A
    pred_M[:, :, 0, 1] = pred_M_raw[:, :, 1]  # B
    pred_M[:, :, 0, 2] = pred_M_raw[:, :, 2]  # C
    pred_M[:, :, 1, 0] = pred_M_raw[:, :, 1]  # B
    pred_M[:, :, 1, 1] = pred_M_raw[:, :, 3]  # D
    pred_M[:, :, 1, 2] = pred_M_raw[:, :, 4]  # E
    pred_M[:, :, 2, 0] = pred_M_raw[:, :, 2]  # C
    pred_M[:, :, 2, 1] = pred_M_raw[:, :, 4]  # E
    pred_M[:, :, 2, 2] = pred_M_raw[:, :, 5]  # F

    faces_counts_indices = np.array(faces_counts).cumsum()
    pred_M_list = torch.tensor_split(pred_M, tuple(faces_counts_indices[:-1]), dim=1)
    return pred_M_list
