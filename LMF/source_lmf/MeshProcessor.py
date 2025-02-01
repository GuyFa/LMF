import json
import time
import warnings
from scipy.sparse import load_npz
import os
import trimesh
from easydict import EasyDict
import numpy as np
import igl
from utils.net_test_info_utils import Timinigs


warnings.filterwarnings("ignore")

NUM_SAMPLES = 1024


class MeshProcessor:
    """
    Extracts all preprocessing-related data (sample points  for pointnet; wave-kernel-signature, etc.)
    """

    def __init__(
        self,
        vertices,
        faces,
        ttype,
        source_dir=None,
    ):
        """
        :param vertices:
        :param faces:
        :param ttype: the torch data type to use (float, half, double)
        :param source_dir: the directory to load the preprocessed data from; if given, will try to load the data before computing, if not given, always compute
        """

        self.ttype = ttype
        self.num_samples = NUM_SAMPLES
        self.vertices = vertices.squeeze()
        self.faces = faces.squeeze()
        self.is_already_cached = True
        self.original_vertices = None
        self.alignment_scale = None
        self.normals = igl.per_vertex_normals(self.vertices, self.faces)
        self.samples = EasyDict()
        self.timings = Timinigs()
        self.samples.xyz = None
        self.samples.normals = None
        self.points_and_normals = None
        self.source_dir = source_dir

    @staticmethod
    def meshprocessor_from_directory(
        source_dir,
        ttype,
    ):
        try:
            vertices = np.load(os.path.join(source_dir, "vertices.npy"))
            faces = np.load(os.path.join(source_dir, "faces.npy"))
        except IOError:
            print(os.path.join(source_dir, "vertices.npy"))
            import traceback

            traceback.print_exc()
            raise IOError("Could not load verices or faces!")
        return MeshProcessor(vertices, faces, ttype, source_dir)

    def get_vertices(self):
        return self.vertices

    def get_faces(self):
        return self.faces

    def get_samples(self):
        if self.samples.xyz is None:
            try:
                self.load_samples()
            except Exception:
                self.is_already_cached = False
                self.compute_samples()
                self.save_samples()
        return self.samples

    def load_samples(self):
        if self.samples.xyz is None:
            self.samples.xyz = np.load(os.path.join(self.source_dir, "samples.npy"))
        if self.samples.normals is None:
            self.samples.normals = np.load(os.path.join(self.source_dir, "samples_normals.npy"))

    def save_samples(self):
        os.makedirs(self.source_dir, exist_ok=True)
        np.save(os.path.join(self.source_dir, "samples.npy"), self.samples.xyz)
        np.save(os.path.join(self.source_dir, "samples_normals.npy"), self.samples.normals)

    def compute_samples(self):
        pt_samples, normals_samples, bary = self.sample_points(self.num_samples)
        self.samples.xyz = pt_samples
        self.samples.normals = normals_samples

    def get_timings(self):
        if self.is_already_cached:
            with open(os.path.join(self.source_dir, "timings.json")) as file:
                self.timings = Timinigs(**json.load(file))
        else:
            self.get_alignment_timing()
            with open(os.path.join(self.source_dir, "timings.json"), "w") as f:
                json.dump(self.timings.model_dump(), f, ensure_ascii=False, indent=4)

    def get_centroids(self):
        if self.points_and_normals is None:
            self.compute_centroids()
        return self.points_and_normals

    def compute_centroids(self):
        start = time.time()
        m = trimesh.Trimesh(vertices=self.vertices, faces=self.faces, process=False)
        self.points_and_normals = np.hstack((np.mean(m.triangles, axis=1), m.face_normals))
        end = time.time()
        self.timings.points_and_normals = end - start
        self.get_samples()

    def get_vertices_before_alignment(self):
        if self.original_vertices is None:
            self.original_vertices = np.load(os.path.join(self.source_dir, "original_vertices.npy"))

    def get_alignment_scale(self):
        if self.alignment_scale is None:
            self.alignment_scale = np.load(os.path.join(self.source_dir, "object_scale.npy"))[0]

    def get_alignment_timing(self):
        self.timings.object_alignment_staff = np.load(os.path.join(self.source_dir, "alignment_time.npy"))[0]

    def get_data(self, key, file_type="npy"):
        if key == "samples":
            return self.get_samples().xyz
        elif key == "samples_normals":
            return self.get_samples().normals
        elif key == "vertices":
            return self.vertices
        elif key == "faces":
            return self.faces
        if file_type == "npy":
            return np.load(os.path.join(self.source_dir, f"{key}.npy"))
        elif file_type == "npz":
            return load_npz(os.path.join(self.source_dir, f"{key}.npz"))
        else:
            raise RuntimeError("wrong file type")

    def sample_points(self, n):
        i = 0
        start = time.time()
        while True:
            if i > 0:
                print("attempt number " + str(i) + " at generating real normal samples")
            i += 1
            bary, found_faces, _ = igl.random_points_on_mesh(n, self.vertices, self.faces)
            vert_ind = self.faces[found_faces - 1]
            point_samples = (
                self.vertices[vert_ind[:, 0]] * bary[:, 0:1]
                + self.vertices[vert_ind[:, 1]] * bary[:, 1:2]
                + self.vertices[vert_ind[:, 2]] * bary[:, 2:3]
            )
            normal_samples = (
                self.normals[vert_ind[:, 0]] * bary[:, 0:1]
                + self.normals[vert_ind[:, 1]] * bary[:, 1:2]
                + self.normals[vert_ind[:, 2]] * bary[:, 2:3]
            )
            if not (np.isnan(normal_samples).any()) and not (np.isinf(normal_samples).any()):
                break

        end = time.time()
        self.timings.all_samples = end - start
        return point_samples, normal_samples, bary
