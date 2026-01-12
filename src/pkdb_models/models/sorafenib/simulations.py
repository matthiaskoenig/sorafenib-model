"""Run sorafenib simulation experiments."""
import shutil
from pathlib import Path
from typing import List

from pymetadata.console import console

from pkdb_models.models import sorafenib
from pkdb_models.models.sorafenib.experiments.misc import *
from pkdb_models.models.sorafenib.experiments.studies import *
from pkdb_models.models.sorafenib.helpers import run_experiments

from sbmlutils import log

from sbmlsim.plot import Figure
Figure.legend_fontsize = 8
Figure.fig_dpi = 300

logger = log.get_logger(__name__)

EXPERIMENTS = \
    {
        'studies': [
            Aboualfa2006,
            Andriamanana2013,
            Awada2005,
            Bins2017,
            Duran2007,
            # Ferrario2016,
            Fucile2015,
            Fukudo2014,
            Hornecker2012,
            Huang2017,
            Huh2021,
            Hussaarts2020,
            Ishii2014,
            Mammatas2020,
            Strumberg2005,
            Zimmerman2012,  # children
        ],
        'misc': [
            DoseDependencyExperiment,
            HepaticImpairmentExperiment,
            RenalImpairmentExperiment,
        ]
    }
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"]

def run_simulation_experiments(
    selected: str = None,
    experiment_classes: List = None,
    output_dir: Path = None
) -> None:
    """Run sorafenib simulation experiments."""

    Figure.fig_dpi = 600
    Figure.legend_fontsize = 10

    # Determine which experiments to run
    if experiment_classes is not None:
        experiments_to_run = experiment_classes
        if output_dir is None:
            output_dir = sorafenib.RESULTS_PATH_SIMULATION / "custom_selection"
    elif selected:
        # Using the 'selected' parameter
        if selected not in EXPERIMENTS:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(EXPERIMENTS.keys())}[/red]"
            )
            console.rule(style="red bold")
            return
        experiments_to_run = EXPERIMENTS[selected]
        if output_dir is None:
            output_dir = sorafenib.RESULTS_PATH_SIMULATION / selected
    else:
        console.print("\n[red bold]Error: No experiments specified![/red bold]")
        console.print("[yellow]Use selected='all' or selected='studies' or provide experiment_classes=[...][/yellow]\n")
        return

    # Run the experiments
    run_experiments(experiment_classes=experiments_to_run, output_dir=output_dir)

    # Collect figures into one folder
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")

if __name__ == "__main__":
    """
    Run experiments
    """
    run_simulation_experiments(selected="studies")
