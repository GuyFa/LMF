import numpy as np


COMPARISON_TOLERANCE = 1e-13


def extract_k(vertices: np.ndarray, faces: np.ndarray, uvs: np.ndarray):
    paramc = uvs[:, 0] + 1j * uvs[:, 1]

    i_1 = faces[:, 0]
    i_2 = faces[:, 1]
    i_3 = faces[:, 2]

    p1 = vertices[i_1, :]
    p2 = vertices[i_2, :]
    p3 = vertices[i_3, :]

    u1 = p3 - p2
    u2 = p1 - p2

    cross_prod = np.sqrt(np.sum((np.cross(u1, u2, axis=1) ** 2), axis=1))

    u1_length = np.sqrt(np.sum((u1**2), axis=1))
    dot = np.sum((u1 * u2), axis=1)

    triangleAreas = cross_prod / 2

    e1 = u1_length
    e3 = -(dot + 1j * cross_prod) / u1_length
    e2 = -e1 - e3

    t1 = -1j * e1
    t2 = -1j * e2
    t3 = -1j * e3

    f1 = paramc[i_1]
    f2 = paramc[i_2]
    f3 = paramc[i_3]

    fzbarAbs = abs((f1 * t1 + f2 * t2 + f3 * t3) / (4 * triangleAreas))
    fzAbs = abs((f1 * np.conj(t1) + f2 * np.conj(t2) + f3 * np.conj(t3)) / (4 * triangleAreas))

    validIndices = fzAbs != 0

    k = np.zeros((faces.shape[0], 1))
    k[fzAbs == 0, 0] = 1
    k[validIndices, 0] = fzbarAbs[validIndices] / fzAbs[validIndices]
    k = k

    sigma_1 = fzAbs + fzbarAbs
    sigma_2 = fzAbs - fzbarAbs

    return (
        k,
        triangleAreas.reshape((-1, 1)),
        sigma_1.reshape((-1, 1)),
        sigma_2.reshape((-1, 1)),
    )


def compute_E_sd(vertices: np.ndarray, faces: np.ndarray, uvs: np.ndarray):
    k, triangles_areas, sigma_1, sigma_2 = extract_k(vertices, faces, uvs)
    tot_area = np.sum(triangles_areas)
    sigma_1_squared = sigma_1**2
    sigma_2_squared = sigma_2**2
    energy_of_each_triangle = 0.5 * (sigma_1_squared + sigma_2_squared + 1 / sigma_1_squared + 1 / sigma_2_squared)
    e_sd = ((energy_of_each_triangle * triangles_areas) / tot_area).sum()
    return e_sd


def compute_sigmas_from_metric_tensors(metric_tensors):
    assert (
        metric_tensors[:, 0, 1] - metric_tensors[:, 1, 0] < COMPARISON_TOLERANCE
    ).all(), "metric tensors must be symmetric matrices!"
    E, F, G = metric_tensors[:, 0, 0], metric_tensors[:, 0, 1], metric_tensors[:, 1, 1]
    sigma_1 = np.sqrt(0.5 * (E + G + np.sqrt(4 * F**2 + (E - G) ** 2)))
    sigma_2 = np.sqrt(0.5 * (E + G - np.sqrt(4 * F**2 + (E - G) ** 2)))
    return sigma_1, sigma_2
