from pathlib import Path

SORAFENIB_PATH = Path(__file__).parent

MODEL_BASE_PATH = SORAFENIB_PATH / "models" / "results" / "models"
MODEL_PATH = MODEL_BASE_PATH / "sorafenib_body_flat.xml"

RESULTS_PATH = SORAFENIB_PATH / "results"
RESULTS_PATH_SIMULATION = RESULTS_PATH / "simulation"
RESULTS_PATH_FIT = RESULTS_PATH / "fit"

H5_PATH = RESULTS_PATH / "pkdb.h5"

# DATA_PATH_BASE = SORAFENIB_PATH.parents[3] / "pkdb_data" / "studies"
DATA_PATH_BASE = SORAFENIB_PATH / "data"
DATA_PATH_SORAFENIB = DATA_PATH_BASE / "sorafenib"
DATA_PATHS = [
    DATA_PATH_SORAFENIB,
]

