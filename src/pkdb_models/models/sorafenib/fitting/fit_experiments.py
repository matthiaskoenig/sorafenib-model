"""Parameter fit problems for sorafenib."""
from sbmlutils.log import get_logger

from pkdb_models.models.sorafenib import SORAFENIB_PATH, DATA_PATHS

from pkdb_models.models.sorafenib.experiments.studies import *

from typing import Dict, List, Type, Union, Callable

from sbmlsim.experiment import ExperimentRunner, SimulationExperiment
from sbmlsim.fit import FitExperiment, FitMapping

from pkdb_models.models.sorafenib.experiments.studies.ferrario2016 import Ferrario2016
from pkdb_models.models.sorafenib.experiments.studies.fukudo2014 import Fukudo2014
from pkdb_models.models.sorafenib.experiments.studies.hussaarts2020 import Hussaarts2020
from pkdb_models.models.sorafenib.experiments.studies.ishii2014 import Ishii2014

logger = get_logger(__name__)


def fit_experiments_for_filter(
    experiment_classes: Union[
        Type[SimulationExperiment], List[Type[SimulationExperiment]]
    ],
    metadata_filter: Callable,
    use_mapping_weights: bool = True,
    default_weight: bool = None,
) -> Dict[str, List[FitExperiment]]:
    """Fit experiments based on MappingMetaData.

    :param experiment_classes: List of SimulationExperiment class definition
    :param metadata_filter:
    """
    # instantiate objects for filtering of fit mappings
    runner = ExperimentRunner(
        experiment_classes=experiment_classes,
        base_path=SORAFENIB_PATH,
        data_path=DATA_PATHS,
    )

    fit_experiments: Dict[str, List[FitExperiment]] = {}

    for k, experiment_name in enumerate(runner.experiments):
        # print(experiment_name)
        experiment_class = experiment_classes[k]
        experiment = runner.experiments[experiment_name]

        # filter mappings by metadata
        mappings = []
        for fm_key, fit_mapping in experiment.fit_mappings.items():
            if metadata_filter(fm_key, fit_mapping):
                mappings.append(fm_key)

        if mappings:
            # add fit experiment from filtered mappings
            fit_experiments[experiment_name] = [
                FitExperiment(
                    experiment=experiment_class,
                    mappings=mappings,
                    weights=default_weight,
                    use_mapping_weights=use_mapping_weights,
                )
            ]

    return fit_experiments


def filter_empty(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Return all experiments/mappings."""
    return True


class Ishi2014:
    pass


experiment_classes = [
    Aboualfa2006,
    Andriamanana2013,
    Awada2005,
    Bins2017,  # TODO: exclude RIF !
    Duran2007,
    Ferrario2016,
    Fucile2015,
    Fukudo2014,
    Hornecker2012,
    Huang2017,
    Huh2021,
    Hussaarts2020,  # TODO: exclude PROB !
    Ishii2014,
    Mammatas2020,
    Strumberg2005,
    Zimmerman2012,
]


logger.info("--- ALL ---")
fitexp_all = fit_experiments_for_filter(
    experiment_classes,
    metadata_filter=filter_empty,
)
logger.info(fitexp_all)
