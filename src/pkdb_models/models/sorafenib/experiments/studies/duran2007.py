from copy import deepcopy
from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models import sorafenib

from pkdb_models.models.sorafenib.experiments.base_experiment import (
    SorafenibSimulationExperiment
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.sorafenib.helpers import run_experiments


class Duran2007(SorafenibSimulationExperiment):
    """Simulation experiment of Duran2007."""
    doses = [200, 400]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)

            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                dset.unit_conversion("mean", 1 / self.Mr.sor)
                dsets[label] = dset

        # print(dsets.keys())
        # print(dsets)

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for dose in self.doses:
            tc_first_day = Timecourse(
                start=0,
                end=24 * 60,  # [min] (24 hours)
                steps=1000,
                changes={
                    **self.default_changes(),
                    "PODOSE_sor": Q_(dose, 'mg'),
                }
            )
            tc_dosing = Timecourse(
                start=0,
                end=24 * 60,  # [min] (24 hours)
                steps=1000,
                changes={
                    "PODOSE_sor": Q_(dose, 'mg'),
                }
            )

            tcsims[f"sor_po_{dose}_multi"] = TimecourseSim(
                [tc_first_day] + [deepcopy(tc_dosing) for _ in range(4)],
                time_offset=-4 * 24 * 60
            )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for dose in [200, 400]:
            mappings[f"fm_sor_po_{dose}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"SOF{dose}_sor",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,  task=f"task_sor_po_{dose}_multi", xid="time", yid="[Cve_sor]",
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        name = "Fig2"
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="hour"), legend=True)
        plots[0].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])

        # simulation
        for dose in self.doses:
            # simulation
            plots[0].add_data(
                task=f"task_sor_po_{dose}_multi",
                xid="time",
                yid="[Cve_sor]",
                label=f"Sim {dose} [mg])",
                color=self.color_for_dose(dose),
            )

            # data
            dataset = f"SOF{dose}_sor"
            plots[0].add_data(
                dataset=dataset,
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
                label=f"{dose} [mg] ({dataset})",
                color=self.color_for_dose(dose),
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Duran2007.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Duran2007, output_dir="Duran2007")
