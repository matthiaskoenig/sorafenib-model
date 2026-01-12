from copy import deepcopy
from typing import Dict

from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.sorafenib.experiments.base_experiment import (
    SorafenibSimulationExperiment,
)
from pkdb_models.models.sorafenib.helpers import run_experiments


class DoseDependencyExperiment(SorafenibSimulationExperiment):
    """Tests multi dose oral sorafenib."""

    doses = [
        100, 200, 400, 800, 1600, 3200
    ]

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for dose in self.doses:
            tc0 = Timecourse(
                start=0,
                end=24*60,
                steps=200,
                changes={
                    **self.default_changes(),
                    "PODOSE_sor": Q_(dose, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24*60,
                steps=100,
                changes={
                    "PODOSE_sor": Q_(dose, "mg"),
                },
            )

            tcsims[f"sor_po{dose}_single"] = TimecourseSim(
                tc0
            )
            tcsims[f"sor_po{dose}_multi"] = TimecourseSim(
                [tc0] + [deepcopy(tc1) for _ in range(10)]
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:

        figs: Dict[str, Figure] = {}
        for suffix in ["single", "multi"]:
            name = f"Oral sorafenib: {suffix.title()} dose"
            fig = Figure(
                experiment=self,
                sid=f"Fig_application_{suffix}",
                num_rows=4,
                num_cols=2,
                name=name,
            )
            plots = fig.create_plots(xaxis=Axis("time", unit="day"), legend=True)

            sids = [
                "[Cve_sor]",
                "[Cve_m2]",
                "[Cve_sg]",
                "Afeces_sor",
                "Aurine_sg",
                "Afeces_sg",
                "Cve_sg_sor",
                "Cve_m2_sor",
            ]

            for k, sid in enumerate(sids):
                plots[k].set_yaxis(
                    self.labels[sid],
                    unit=self.units[sid],
                )

                for dose in self.doses:
                    task_id = f"task_sor_po{dose}_{suffix}"
                    # simulation
                    plots[k].add_data(
                        task=task_id,
                        xid="time",
                        yid=sid,
                        label=f"SOR {dose} mg PO",
                        color=self.color_for_dose(dose)
                    )

            figs[fig.sid] = fig

        return figs


if __name__ == "__main__":
    run_experiments(DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__)
