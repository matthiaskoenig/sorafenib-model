"""Sorafenib parameter fitting."""
import logging
from pathlib import Path

import itertools
from typing import List, Dict, Tuple

from sbmlsim.fit import FitParameter, FitExperiment
from sbmlsim.fit.result import OptimizationResult
from sbmlsim.fit.optimization import OptimizationProblem
from sbmlsim.fit.analysis import OptimizationAnalysis
from sbmlsim.fit.runner import run_optimization
from sbmlsim.fit.options import *
from sbmlsim.fit.sampling import SamplingType

from pkdb_models.models.sorafenib.fitting.fit_experiments import (
    fitexp_all,
)
from pkdb_models.models.sorafenib.fitting.parameters import (
    parameters_all,
)

from pkdb_models.models.sorafenib import (
    RESULTS_PATH_FIT,
    SORAFENIB_PATH,
    DATA_PATHS,
)

logger = logging.getLogger(__name__)

fit_kwargs = {
    # optimization settings
    "residual": ResidualType.NORMALIZED,
    "loss_function": LossFunctionType.LINEAR,
    "weighting_curves": [
        WeightingCurvesType.MAPPING,  # user defined weights # FIXME
        WeightingCurvesType.POINTS,  # number of points
    ],
    "weighting_points": WeightingPointsType.ERROR_WEIGHTING,  # no errors weighed with CV=0.5
    # additional integrator settings
    "variable_step_size": True,
    "relative_tolerance": 1e-6,
    "absolute_tolerance": 1e-6,
    # serial optimization
    # "serial": False,
}


def create_optimization_problem(
    fit_experiments: List[FitExperiment], opid: str, parameters: List[FitParameter]
) -> OptimizationProblem:
    op = OptimizationProblem(
        opid=opid,
        fit_experiments=fit_experiments,
        fit_parameters=parameters,
        base_path=SORAFENIB_PATH,
        data_path=DATA_PATHS,
    )
    return op


def fitlsq(op, **kwargs) -> Tuple[OptimizationResult, OptimizationProblem]:
    """Local least square fitting."""
    opt_res = run_optimization(
        problem=op,
        seed=1238,
        algorithm=OptimizationAlgorithmType.LEAST_SQUARE,
        # parameters for least square optimization
        sampling=SamplingType.LOGUNIFORM_LHS,
        diff_step=0.05,
        # diff_step=0.05,
        # ftol=1e-10,
        # xtol=1e-10,
        # gtol=1e-10,
        **kwargs
    )
    return opt_res, op


def fitde(op, **kwargs) -> Tuple[OptimizationResult, OptimizationProblem]:
    """Global differential evolution fitting."""
    opt_res = run_optimization(
        problem=op,
        seed=1234,
        algorithm=OptimizationAlgorithmType.DIFFERENTIAL_EVOLUTION,
        **kwargs
    )
    return opt_res, op


class OptimizationStrategy(Enum):
    """Strategy for fitting.

    Either fit all experiments together or individually.
    """

    ALL = 1  # fit all experiments together, i.e, one parameter set
    SINGLE = (
        2  # fit individual experiments, i.e. set of parameters for every experiment
    )


class FitMethod(Enum):
    """Method for fitting."""

    LSQ = 1
    DE = 2


class FitExperimentSubset(Enum):
    """Subset of fit experiments for fitting."""

    ALL = 1


def fit_sorafenib(
    optimization_strategy: OptimizationStrategy,
    fit_method: FitMethod,
    fit_experiments: List[FitExperiment],
    parameters: List[FitParameter],
) -> Dict[str, Tuple[OptimizationResult, OptimizationProblem]]:

    if not isinstance(optimization_strategy, OptimizationStrategy):
        raise ValueError
    if not isinstance(fit_method, FitMethod):
        raise ValueError

    def fit_op(
        op: OptimizationProblem,
    ) -> Tuple[OptimizationResult, OptimizationProblem]:
        """Wrapper for optimization function."""

        # run optimization
        opt_result: OptimizationResult
        op: OptimizationProblem
        if fit_method == FitMethod.LSQ:
            opt_result, op = fitlsq(op, size=10, n_cores=10, **fit_kwargs)
        elif fit_method == FitMethod.DE:
            opt_result, op = fitde(op, size=10, n_cores=10, **fit_kwargs)

        return opt_result, op

    # store optimization results
    results = {}

    if optimization_strategy == OptimizationStrategy.SINGLE:
        # fit all experiments individually
        for fit_exp in fit_experiments:
            opid = fit_exp.experiment_class.__name__

            op = create_optimization_problem(
                fit_experiments=[fit_exp], opid=opid, parameters=parameters
            )
            results[opid] = fit_op(op=op)

    elif optimization_strategy == OptimizationStrategy.ALL:
        # fit all experiments together
        opid = "all"
        op = create_optimization_problem(
            fit_experiments=fit_experiments, opid=opid, parameters=parameters
        )
        results[opid] = fit_op(op)

    return results


def get_fit_experiments(fit_subset: FitExperimentSubset, study_ids: List[str] = None):
    """Creates a subset of fit experiments from given information."""
    if not isinstance(fit_subset, FitExperimentSubset):
        raise ValueError

    if fit_subset == FitExperimentSubset.ALL:
        fitexp_dict = fitexp_all

    if study_ids:
        fit_experiments = [fitexp_dict[sid] for sid in study_ids]
    else:
        fit_experiments = [exp for exp in fitexp_dict.values()]

    # reduce list of lists
    fit_experiments = list(itertools.chain(*fit_experiments))
    return fit_experiments


def get_fit_parameters(fit_subset: FitExperimentSubset) -> List[FitParameter]:
    """Creates a subset of fit experiments from given information."""
    if not isinstance(fit_subset, FitExperimentSubset):
        raise ValueError

    if fit_subset == FitExperimentSubset.ALL:
        parameters = parameters_all

    return parameters


def run_fit_subset(
    fit_experiments: List[FitExperiment],
    parameters: List[FitParameter],
    fit_method: FitMethod,
    output_name: str,
    output_dir: Path,
):
    """Fits subset of data.

    Used for iterative fitting of parameters.
    """
    results_all: Dict[str, Tuple[OptimizationResult, OptimizationProblem]] = fit_sorafenib(
        fit_experiments=fit_experiments,
        parameters=parameters,
        optimization_strategy=OptimizationStrategy.ALL,
        fit_method=fit_method,
    )

    # parameters for plots
    mpl_parameters = {
        # 'axes.labelsize': 12,
        # 'axes.labelweight': "bold",
    }
    for key, (opt_result, op) in results_all.items():
        # create figures and outputs
        opt_analysis = OptimizationAnalysis(
            opt_result=opt_result,
            op=op,
            output_name=output_name,
            output_dir=output_dir,
            show_plots=True,
            show_titles=False,
            **fit_kwargs
        )
        opt_analysis.run(mpl_parameters=mpl_parameters)

    return results_all


if __name__ == "__main__":
    """
    Fitting of sorafenib parameters.
    """
    fit_definition = [
        [
            FitExperimentSubset.ALL,
            FitMethod.LSQ,
            # FitMethod.DE,
            "lsq_all",
        ],
    ]

    results: Dict[str, Tuple[OptimizationResult, OptimizationProblem]]
    results = {}
    for (fit_subset, fit_method, output_name) in fit_definition:

        fit_experiments = get_fit_experiments(fit_subset=fit_subset)
        parameters = get_fit_parameters(fit_subset=fit_subset)

        results = run_fit_subset(
            fit_experiments=fit_experiments,
            parameters=parameters,
            fit_method=fit_method,
            output_name=output_name,
            output_dir=RESULTS_PATH_FIT,
        )
