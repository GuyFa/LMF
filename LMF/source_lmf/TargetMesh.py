import torch
from easydict import EasyDict

import MeshProcessor


class TargetMesh:
    """
    datastructure for the target mappings, in batched form (one source can be mapped to multiple targets)
    """

    def __init__(self, target_ind, target_dir, extra_target_fields, ttype):
        self.target_dir = target_dir
        self.target_ind = target_ind

        self.target_tensors = None

        self.__extra_keys = extra_target_fields
        self.__loaded_data = {}
        self.ttype = ttype
        self.mesh_processor = None

    def get_loaded_data(self, key: str):
        return self.__loaded_data.get(key)

    def to(self, device):
        for attr in self.target_tensors.__dict__.keys():
            t = getattr(self.target_tensors, attr)
            setattr(self.target_tensors, attr, t.to(device))
        for key in self.__loaded_data.keys():
            self.__loaded_data[key] = self.__loaded_data[key].to(device)
        return self

    def load(self):
        self.target_tensors = EasyDict()
        self.mesh_processor = MeshProcessor.MeshProcessor.meshprocessor_from_directory(
            self.target_dir,
            self.ttype,
        )
        self.target_tensors.vertices = self.mesh_processor.get_data("vertices")

        for attr in self.target_tensors.__dict__.keys():
            t = getattr(self.target_tensors, attr)
            setattr(self.target_tensors, attr, torch.from_numpy(t).type(self.ttype))

    def get_vertices(self):
        return self.target_tensors.vertices
