import numpy as np
import igl
import torch


def weighted_mse_loss(input, target, weights, norm_axis, percentage=1, order="increase"):
    # the weights can be unnormalized - we normalize it here
    diff = target - input
    errors = torch.linalg.norm(diff, axis=norm_axis) ** 2
    if percentage == 1:
        loss = torch.sum(weights * errors) / weights.sum()
        return loss
    else:
        sort_indices = np.argsort(errors)
        percentage_count = np.round(percentage * errors.shape[0])
        if order == "increase":
            indices = sort_indices[:percentage_count]
        else:
            indices = sort_indices[percentage_count:]
        loss = torch.sum(weights[indices] * errors[indices]) / weights[indices].sum()
        return loss


def compute_matrix_loss(source, gt_matrix_list, pred_matrix_list):
    matrix_loss = torch.tensor(0.0, device="cuda:0")
    for i in range(len(source)):
        source_areas = torch.from_numpy(
            igl.doublearea(source[i].mesh_processor.vertices, source[i].mesh_processor.faces)
        ).to("cuda:0")

        cur_matrix_loss = weighted_mse_loss(
            gt_matrix_list[i].squeeze(),
            pred_matrix_list[i].squeeze(),
            source_areas,
            (1, 2),
        )
        matrix_loss += (1 / len(source)) * cur_matrix_loss
    return matrix_loss


def log_mat(matrices_tensor):
    # Assuming PSD matrices
    if matrices_tensor.shape[0] > 1.5e6 and matrices_tensor.shape[0] < 3e6:
        eig_vals_1, eig_vecs_1 = torch.linalg.eigh(matrices_tensor[:1500000, :, :])
        eig_vals_2, eig_vecs_2 = torch.linalg.eigh(matrices_tensor[1500000:, :, :])
        eig_vals = torch.vstack((eig_vals_1, eig_vals_2))
        eig_vecs = torch.vstack((eig_vecs_1, eig_vecs_2))
    elif matrices_tensor.shape[0] > 3e6 and matrices_tensor.shape[0] < 4.5e6:
        eig_vals_1, eig_vecs_1 = torch.linalg.eigh(matrices_tensor[:1500000, :, :])
        eig_vals_2, eig_vecs_2 = torch.linalg.eigh(matrices_tensor[1500000:3000000, :, :])
        eig_vals_3, eig_vecs_3 = torch.linalg.eigh(matrices_tensor[3000000:, :, :])
        eig_vals = torch.vstack((eig_vals_1, eig_vals_2, eig_vals_3))
        eig_vecs = torch.vstack((eig_vecs_1, eig_vecs_2, eig_vecs_3))
    else:
        eig_vals, eig_vecs = torch.linalg.eigh(matrices_tensor)
    log_eig_vals_mat = torch.zeros_like(matrices_tensor)
    log_eig_vals_mat[:, 0, 0] = eig_vals[:, 0].log()
    log_eig_vals_mat[:, 1, 1] = eig_vals[:, 1].log()
    log_eig_vals_mat[:, 2, 2] = eig_vals[:, 2].log()
    res = torch.einsum(
        "abc,acd->abd",
        torch.einsum("abc,acd->abd", torch.real(eig_vecs), log_eig_vals_mat),
        torch.real(eig_vecs).transpose(1, 2),
    )

    return res


def exp_mat(matrices_tensor: torch.Tensor):
    if matrices_tensor.shape[0] > 1.5e6 and matrices_tensor.shape[0] < 3e6:
        eig_vals_1, eig_vecs_1 = torch.linalg.eigh(matrices_tensor[:1500000, :, :])
        eig_vals_2, eig_vecs_2 = torch.linalg.eigh(matrices_tensor[1500000:, :, :])
        eig_vals = torch.vstack((eig_vals_1, eig_vals_2))
        eig_vecs = torch.vstack((eig_vecs_1, eig_vecs_2))
    elif matrices_tensor.shape[0] > 3e6 and matrices_tensor.shape[0] < 4.5e6:
        eig_vals_1, eig_vecs_1 = torch.linalg.eigh(matrices_tensor[:1500000, :, :])
        eig_vals_2, eig_vecs_2 = torch.linalg.eigh(matrices_tensor[1500000:3000000, :, :])
        eig_vals_3, eig_vecs_3 = torch.linalg.eigh(matrices_tensor[3000000:, :, :])
        eig_vals = torch.vstack((eig_vals_1, eig_vals_2, eig_vals_3))
        eig_vecs = torch.vstack((eig_vecs_1, eig_vecs_2, eig_vecs_3))
    else:
        eig_vals, eig_vecs = torch.linalg.eigh(matrices_tensor)
    exp_eig_vals_mat = torch.zeros_like(matrices_tensor)
    exp_eig_vals_mat[:, 0, 0] = eig_vals[:, 0].exp()
    exp_eig_vals_mat[:, 1, 1] = eig_vals[:, 1].exp()
    exp_eig_vals_mat[:, 2, 2] = eig_vals[:, 2].exp()
    res = torch.einsum(
        "abc,acd->abd",
        torch.einsum("abc,acd->abd", torch.real(eig_vecs), exp_eig_vals_mat),
        torch.real(eig_vecs).transpose(1, 2),
    )

    return res
