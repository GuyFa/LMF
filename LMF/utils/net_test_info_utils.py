from typing import List, Optional
import numpy as np
import torch
from pydantic import BaseModel, ConfigDict, computed_field


class Timinigs(BaseModel):
    object_alignment_staff: float = 0
    all_samples: float = 0
    points_and_normals: float = 0
    move_data_to_gpu: float = 0
    predict_mt: float = 0
    compute_exp_mt: float = 0
    run_GIF: float = 0

    @computed_field
    @property
    def total_time(self) -> float:
        return (
            self.object_alignment_staff
            + self.all_samples
            + self.points_and_normals
            + self.move_data_to_gpu
            + self.predict_mt
            + self.compute_exp_mt
            + self.run_GIF
        )


class TestStep3DMetricTensorOutput(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)
    target_V: List[np.ndarray]
    source_V: List[np.ndarray]
    total_loss: torch.Tensor
    T: List[np.ndarray]
    source_ind: List[str]
    target_inds: List[str]
    batch_size: int
    pred_V: List[np.ndarray]
    pred_flips_before_fix: np.ndarray
    pred_flips_after_fix: np.ndarray
    IDT_stats: List[np.ndarray]
    timer: Timinigs
    edge_compatibility_data: Optional[tuple] = None
