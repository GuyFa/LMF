import random
from typing import Literal

from scipy.linalg import orthogonal_procrustes
from torch.utils.data import Dataset
import numpy as np
import os
from TargetMesh import TargetMesh
from SourceMesh import SourceMesh
import time
from os.path import join
from pathlib import Path
import logging
import json


logging.basicConfig(level=logging.DEBUG, filename="./exception.log")


class DeformationDataset(Dataset):
    """
    Main dataset. Each sample in it is a source <---> target pair. Note this dataset return a custom Batch class object instead of a
    tensor, which is already batched.
    """

    def num_meshes(self):
        return len(self.file_names)

    def num_pairs(self):
        return self.len_pairs

    def __init__(
        self,
        path,
        data_file,
        source_keys,
        target_keys,
        ttype,
        args,
        dataset_type: Literal["train", "validation", "test"],
    ):
        """
        :param s_and_t: list of tuples, each tuple is two indices into the file_names, giving a source/target pair
        :param max_source_batch: max batch_size. Note batches cannot have more than one unique source mesh.
        :param batch_class: which class to use to create batch objects
        """
        SHUFFLE_TARGETS_OF_SINGLE_SOURCE = True
        with open(os.path.join(path, data_file)) as file:
            data = json.load(file)
            s_and_t = data["pairs"]
        self.ttype = ttype
        self.dataset_type = dataset_type
        self.unique_source = False

        self.len_pairs = len(s_and_t)
        if SHUFFLE_TARGETS_OF_SINGLE_SOURCE:
            random.shuffle(s_and_t)

        self.source_and_target = None
        source_target = {}
        for val in s_and_t:
            s = val[0]
            t = val[1]
            if s not in source_target:
                source_target[s] = t
        if len(source_target) == 1:
            # This flag will avoid reloading the source at each iteration, since the source is always the same.
            self.unique_source = True
        self.source_and_target = []
        for s, ts in source_target.items():
            self.source_and_target.append((s, ts))

        self.source_keys = source_keys
        self.s_and_t = s_and_t
        self.target_keys = target_keys
        self.args = args
        self.directory = path
        self.source = None
        self.target = None

    def __len__(self):
        if self.source_and_target is None:
            return self.len_pairs
        return len(self.source_and_target)

    def get_item_default(self, ind, verbose=False):
        # Single source single target
        source_index, target_index = self.source_and_target[ind]

        assert Path(source_index).suffix in [".obj", ".off", ".ply"] and Path(target_index).suffix in [
            ".obj",
            ".off",
            ".ply",
        ]

        target_directory_name = Path(join(self.directory, target_index)).as_posix()[:-4]
        if not os.path.exists(join(target_directory_name, "vertices.npy")) and not os.path.exists(
            join(target_directory_name, "faces.npy")
        ):
            source_scale = self.obj_to_npy(Path(join(self.directory, source_index)), True)
            self.obj_to_npy(Path(join(self.directory, target_index)), False, source_scale)

        source_index = source_index[:-4]
        target_index = target_index[:-4]

        # ==================================================================
        # LOAD SOURCE
        if self.source is not None and self.unique_source:
            source = self.source
        else:
            source = SourceMesh(
                source_index,
                join(self.directory, source_index),
                self.source_keys,
                self.ttype,
                target_dir=join(self.directory, target_index),
            )
            source.load()
            self.source = source

        # ==================================================================
        # LOAD TARGET
        target = TargetMesh(
            target_index,
            join(self.directory, target_index),
            self.target_keys,
            self.ttype,
        )
        target.load()
        return source, target

    def obj_to_npy(self, path, is_source=True, scale=None):
        assert is_source or scale is not None
        import igl

        directory_name = path.as_posix()[:-4]
        if not os.path.exists(join(directory_name, "vertices.npy")) and not os.path.exists(
            join(directory_name, "faces.npy")
        ):
            os.makedirs(directory_name, exist_ok=True)
            mesh = igl.read_triangle_mesh(path.as_posix())
            vertices, faces = mesh[0].copy(), mesh[1].copy()
            np.save(join(directory_name, "original_vertices.npy"), vertices)
            if is_source:
                bnd = igl.boundary_loop(faces)  # Assuming boundary vertices are known before the run
                start_time = time.time()
                scale = 1 / np.sqrt(np.sum(igl.doublearea(vertices, faces)))
                vertices = scale * np.array(vertices.tolist())
                if not self.args.no_rotational_alignment:
                    vertices = rotation_boundary_PCA(vertices, bnd)
                source_mesh_centroid = (vertices.max(axis=0) + vertices.min(axis=0)) / 2
                vertices -= source_mesh_centroid
                end_time = time.time()
                np.save(join(directory_name, "object_scale.npy"), np.array([scale]))
                np.save(
                    join(directory_name, "alignment_time.npy"),
                    np.array([end_time - start_time]),
                )
            if not is_source:
                vertices = scale * np.array(vertices.tolist())
                vertices -= vertices.mean(axis=0)
            np.save(join(directory_name, "vertices.npy"), vertices)
            np.save(join(directory_name, "faces.npy"), faces)
        return scale

    def __getitem__(self, ind, verbose=False):
        start = time.time()
        data_sample = self.get_item_default(ind)
        if verbose:
            print(f"DATALOADER : loaded sample in {time.time() - start}")

        return data_sample

    def get_point_dim(self):
        return self[0][0].get_point_dim()


def PCA(X, num_components):
    # Step-1
    X_meaned = X - np.mean(X, axis=0)

    # Step-2
    cov_mat = np.cov(X_meaned, rowvar=False)

    # Step-3
    eigen_values, eigen_vectors = np.linalg.eigh(cov_mat)

    # Step-4
    sorted_index = np.argsort(eigen_values)[::-1]
    # sorted_eigenvalue = eigen_values[sorted_index]
    sorted_eigenvectors = eigen_vectors[:, sorted_index]

    # Step-5
    eigenvector_subset = sorted_eigenvectors[:, 0:num_components]

    # Step-6
    X_reduced = np.dot(eigenvector_subset.transpose(), X_meaned.transpose()).transpose()

    return X_reduced


def rotation_boundary_PCA(V, bnd):
    bnd_PCA_V = PCA(V[bnd, :], 3)

    R, _ = orthogonal_procrustes(V[bnd, :], bnd_PCA_V)

    V_local = V @ R

    return V_local
