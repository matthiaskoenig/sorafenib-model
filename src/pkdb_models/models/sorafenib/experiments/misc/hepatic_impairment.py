from copy import deepcopy
from typing import Dict

from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.sorafenib.experiments.base_experiment import (
    SorafenibSimulationExperiment,
)
from pkdb_models.models.sorafenib.helpers import run_experiments


class HepaticImpairmentExperiment(SorafenibSimulationExperiment):
    """Tests hepatic impairment oral sorafenib."""

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for cirrhosis_key in self.cirrhosis_map:

            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=200,
                changes={
                    **self.default_changes(),
                    "PODOSE_sor": Q_(200, "mg"),
                    "f_cirrhosis": Q_(
                        self.cirrhosis_map[cirrhosis_key], "dimensionless"
                    ),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=200,
                changes={
                    "PODOSE_sor": Q_(200, "mg"),
                },
            )

            tcsims[f"sor_po_{cirrhosis_key}"] = TimecourseSim(
                [tc0] + [deepcopy(tc1) for _ in range(10)],
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:

        # name = "Hepatic impairment (cirrhosis)"
        name = None
        fig = Figure(
            experiment=self,
            sid=f"Fig_application",
            num_rows=2,
            num_cols=3,
            name=name,
        )
        plots = fig.create_plots(xaxis=Axis("time", unit="day"), legend=True)

        sids = [
            "[Cve_sor]",
            "[Cve_m2]",
            "[Cve_sg]",
            "Afeces_sor",
            "Afeces_sg",
            "Aurine_sg",
            # "Cve_sg_sor",
            # "Cve_m2_sor",
        ]

        for k, sid in enumerate(sids):
            self.plot_grouped_species(sid, plots[k])

        return {
            fig.sid: fig,
        }

    def plot_grouped_species(self, sid: str, plot: Plot):

        plot.set_yaxis(
            self.labels[sid],
            unit=self.units[sid],
        )

        for cirrhosis_key in self.cirrhosis_map:
            task_id = f"task_sor_po_{cirrhosis_key}"
            # simulation
            plot.add_data(
                task=task_id,
                xid="time",
                yid=sid,
                label=f"{cirrhosis_key}",
                color=self.cirrhosis_colors[cirrhosis_key]
            )


if __name__ == "__main__":
    run_experiments(HepaticImpairmentExperiment, output_dir=HepaticImpairmentExperiment.__name__)
