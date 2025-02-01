from typing import List

import cupy
from pydantic import BaseModel
import numpy as np


def release_memory():
    mempool = cupy.get_default_memory_pool()
    pinned_mempool = cupy.get_default_pinned_memory_pool()
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()


class LossManager(BaseModel):
    data_type: str = ""
    total_loss: List[float] = []

    def print_loss(self):
        print("\n" + self.data_type + " total loss track:\n" + str(self.total_loss))
        return

    def plot_losses(self, ax):
        ax.plot(
            np.arange(len(self.total_loss)),
            np.array(self.total_loss),
            label=self.data_type + " total loss",
        )
