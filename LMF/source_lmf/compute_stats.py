import numpy as np
import numpy.random

from utils.func_utils import extract_k


def compute_stats_from_V(vertices, faces, uvs):
    k_per_triangle, triangles_areas, sigma_1, sigma_2 = extract_k(vertices, faces, uvs)
    sigma_1_squared = sigma_1**2
    sigma_2_squared = sigma_2**2
    e_sd_per_triangle = 0.5 * (sigma_1_squared + sigma_2_squared + 1 / sigma_1_squared + 1 / sigma_2_squared)

    e_sd = np.sum(triangles_areas * e_sd_per_triangle) / np.sum(triangles_areas)
    k_energy = np.sum(triangles_areas * k_per_triangle) / np.sum(triangles_areas)

    flips = (sigma_2 <= 0).sum()
    return (
        e_sd,
        k_energy,
        flips,
    )
