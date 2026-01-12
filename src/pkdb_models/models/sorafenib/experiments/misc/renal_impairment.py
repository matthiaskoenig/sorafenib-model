from copy import deepcopy
from typing import Dict

from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.sorafenib.experiments.base_experiment import (
    SorafenibSimulationExperiment,
)
from pkdb_models.models.sorafenib.helpers import run_experiments


class RenalImpairmentExperiment(SorafenibSimulationExperiment):
    """Tests renal impairment oral sorafenib."""

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for renal_key in self.renal_map:

            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=200,
                changes={
                    **self.default_changes(),
                    "PODOSE_sor": Q_(200, "mg"),
                    "KI__f_renal_function": Q_(
                        self.renal_map[renal_key], "dimensionless"
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

            tcsims[f"sor_po_{renal_key}"] = TimecourseSim(
                [tc0] + [deepcopy(tc1) for _ in range(10)],
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:

        # name = "Renal impairment"
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

        for renal_key in self.renal_map:
            task_id = f"task_sor_po_{renal_key}"
            # simulation
            plot.add_data(
                task=task_id,
                xid="time",
                yid=sid,
                label=f"{renal_key}",
                color=self.renal_colors[renal_key]
            )


if __name__ == "__main__":
    run_experiments(RenalImpairmentExperiment, output_dir=RenalImpairmentExperiment.__name__)
