import time

import torch

import MeshProcessor


class SourceMesh:
    """
    datastructure for the source mesh to be mapped
    """

    def __init__(
        self,
        source_ind,
        source_dir,
        extra_source_fields,
        ttype,
        target_dir=None,
    ):
        self.source_ind = source_ind
        self.source_dir = source_dir
        self.centroids_and_normals = None
        self.__extra_keys = extra_source_fields
        self.__loaded_data = {}
        self.__ttype = ttype
        self.mesh_processor = None
        self.target_dir = target_dir

    def get_vertices(self):
        return self.source_vertices

    def get_loaded_data(self, key: str):
        return self.__loaded_data.get(key)

    def get_source_triangles(self):
        return self.mesh_processor.get_faces()

    def to(self, device):
        start = time.time()

        self.centroids_and_normals = self.centroids_and_normals.to(device)
        for key in self.__loaded_data.keys():
            self.__loaded_data[key] = self.__loaded_data[key].to(device)

        end = time.time()
        self.mesh_processor.timings.move_data_to_gpu += end - start

        return self

    def __init_from_mesh_data(self):
        assert self.mesh_processor is not None
        self.source_vertices = torch.from_numpy(self.mesh_processor.get_vertices()).type(self.__ttype)

        # Load input to MLP
        self.mesh_processor.get_centroids()

        self.mesh_processor.get_timings()

        start = time.time()

        centroid_points_and_normals = self.mesh_processor.points_and_normals
        self.centroids_and_normals = torch.from_numpy(centroid_points_and_normals).type(self.__ttype)

        # Essentially here we load pointnet data and apply the same preprocessing
        for key in self.__extra_keys:
            data = self.mesh_processor.get_data(key)
            data = torch.from_numpy(data)
            data = data.unsqueeze(0).type(self.__ttype)

            self.__loaded_data[key] = data

        end = time.time()
        self.mesh_processor.timings.move_data_to_gpu = end - start

    def load(self):
        self.mesh_processor = MeshProcessor.MeshProcessor.meshprocessor_from_directory(
            self.source_dir,
            self.__ttype,
        )
        self.__init_from_mesh_data()

    def get_point_dim(self):
        return self.centroids_and_normals.shape[1]

    def get_centroids_and_normals(self):
        return self.centroids_and_normals
