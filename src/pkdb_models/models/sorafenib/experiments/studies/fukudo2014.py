from copy import deepcopy
from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models import sorafenib

from pkdb_models.models.sorafenib.experiments.base_experiment import (
    SorafenibSimulationExperiment
)
from pkdb_models.models.sorafenib.helpers import run_experiments


class Fukudo2014(SorafenibSimulationExperiment):
    """Simulation experiment of Fukudo2014"""
    doses = [800]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                dset.unit_conversion("mean", 1 / self.Mr.sor)
                dsets[f"sor_{label}"] = dset

        # print(dsets.keys())
        # print(dsets)

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        # multiple dosing
        #dose = 800  # [mg] twice daily
        Q_ = self.Q_
        tcsims = {}
        for dose in self.doses:
          tc0 = Timecourse(
            start=0,
            end=12 * 60,  # [min]
            steps=200,
            changes={
                **self.default_changes(),
                "PODOSE_sor": Q_(dose/2.0, "mg"),
            },
        )
        tc1 = Timecourse(
            start=0,
            end=12 * 60,  # [min]
            steps=200,
            changes={
                "PODOSE_sor": Q_(dose/2.0, "mg"),
            },
        )

        tcsims[f"sor_po800"] = TimecourseSim(
            [tc0] + [deepcopy(tc1) for _ in range(14)],  # 7 days pretreatment
            time_offset=-7 * 24 * 60
        )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for dose in self.doses:
          mappings[f"fm_sor_po{dose}"] = FitMapping(
            self,
            reference=FitData(
                self,
                dataset=f"sor_RCC_{dose}",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
            ),
            observable=FitData(
                self, task=f"task_sor_po800", xid="time", yid="[Cve_sor]",
            ),
        )
        return mappings


    def figures(self) -> Dict[str, Figure]:
        name = "Fig3"
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=2,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(legend=True)

        # Set x-axis and y-axis labels for each plot
        plots[0].set_xaxis(self.label_time, unit="hr")
        plots[1].set_xaxis(self.label_time, unit="hr", min=-12)
        plots[0].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])
        plots[1].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])

        # simulations
        for dose in self.doses:
          for k in [0, 1]:
            plots[k].add_data(
                task="task_sor_po800",
                xid="time",
                yid="[Cve_sor]",
                label="Sim 800 mg",
                color=self.color_for_dose(dose),
            )

        # Data
        for k in [0, 1]:

            plots[k].add_data(
                dataset=f"sor_RCC_{dose}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"RCC 800 [mg]",
                color=self.color_for_dose(dose),
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Fukudo2014.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Fukudo2014, output_dir="Fukudo2014")
