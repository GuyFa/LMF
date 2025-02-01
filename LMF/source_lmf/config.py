from pydantic import BaseModel


class StatsNetResultsDFCols(BaseModel):
    MODEL_PATH: str = "Model path"
    FACES_COUNT: str = "#Faces"
    VERTICES_COUNT: str = "#Vertices"
    BOUNDARY_VERTICES_COUNT: str = "#Boundary Vertices"
    # GT's stats
    GT_E_SD: str = "GT E_sd"
    GT_K: str = "GT k"
    GT_FLIPS: str = "GT flips"
    # Ours stats
    OURS_E_SD: str = "Ours E_sd"
    OURS_K: str = "Ours k"
    OURS_AFTER_MVC_FIX_FLIPS: str = "Ours #Flips"
    OURS_MEAN_VALUE_FIX_WAS_APPLIED: str = "Ours - Mean Value Fix was Applied"
    OURS_TOTAL_LOSS: str = "Ours Total loss"
    OURS_COMPATIBILITY_ABS_ERROR_MEAN: str = "Ours Compatibility Absolute Error Mean"
    OURS_COMPATIBILITY_ABS_ERROR_STD: str = "Ours Compatibility Absolute Error STD"
    OURS_COMPATIBILITY_REL_ERROR_MEAN: str = "Ours Compatibility Relative Error Mean"
    OURS_COMPATIBILITY_REL_ERROR_STD: str = "Ours Compatibility Relative Error STD"
    OURS_IDT_EDGE_FLIPS_COUNT: str = "Ours IDT Edge Flips Count"
    OURS_IDT_TIME_CONSUMPTION: str = "Ours IDT Time Consumption"
    OURS_NETWORK_TIME_CONSUMPTION: str = "Ours Network Time Consumption"
    OURS_GIF_TIME_CONSUMPTION: str = "Ours GIF Time Consumption"
    OURS_TOTAL_TIME_CONSUMPTION: str = "Ours Total Time Consumption"
    # GIF stats
    GIF_E_SD: str = "GIF E_sd"
    GIF_K: str = "GIF k"
    GIF_AFTER_MVC_FIX_FLIPS: str = "GIF #Flips"
    GIF_MEAN_VALUE_FIX_WAS_APPLIED: str = "GIF - Mean Value Fix was Applied"
    GIF_TIME_CONSUMPTION: str = "GIF Time Consumption"
    # Original CM stats
    CM_E_SD: str = "CM E_SD"
    CM_K: str = "CM k"
    CM_FLIPS: str = "CM #Flips"
    CM_TIME: str = "CM Total Time Consumption"
    CM_ITERS: str = "CM Iterations"
    # CM with ours initialization - stats
    CM_E_SD_WITH_OURS_INIT: str = "CM E_SD with Ours Init"
    CM_TIME_WITH_OURS_INIT: str = "CM Time with Ours Init"
    CM_ITERS_WITH_OURS_INIT: str = "CM Iterations with Ours Init"
    # CM stopped when reaches our energy - stats
    CM_TIME_TO_REACH_OURS_ENERGY: str = "CM time to reach ours E_SD"
    CM_ITERS_TO_REACH_OURS_ENERGY: str = "CM iters to reach ours E_SD"
    CM_E_SD_AFTER_REACHING_OURS_ENERGY: str = "CM E_SD after reaching ours E_SD"
    # CM stopped when reaches our runtime - stats
    CM_ITERS_TO_REACH_OURS_RUNTIME: str = "CM iters to reach ours runtime"
    CM_E_SD_AFTER_REACHING_OURS_RUNTIME: str = "CM E_SD after reaching ours runtime"


def get_cols_names_keys_list(config_cols):
    return list(config_cols.model_dump().keys())


def get_cols_names_values_list(config_cols):
    return list(config_cols.model_dump().values())


class CONFIG:
    def __init__(self):
        self.stats_net_results_df_cols = StatsNetResultsDFCols()


def get_config():
    return CONFIG()
