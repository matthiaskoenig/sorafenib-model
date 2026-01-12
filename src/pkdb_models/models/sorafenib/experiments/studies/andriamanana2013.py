from copy import deepcopy
from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData

from pkdb_models.models.sorafenib.experiments.base_experiment import (
    SorafenibSimulationExperiment
)
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim


from pkdb_models.models.sorafenib.helpers import run_experiments


class Andriamanana2013(SorafenibSimulationExperiment):
    doses = [400]
    """Simulation experiment of Andriamanana2013.
    Multiple dosing with 400 mg BID sorafenib.
    """
    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        fig_id = "Fig4"
        df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)

        for label in df.label.unique():
            df_label = df[df.label == label]
            dset = DataSet.from_df(df_label.copy(), self.ureg)
            dset.unit_conversion("mean", 1 / self.Mr.sor)
            dsets[label] = dset

        # print(dsets.keys())

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for dose in self.doses:
          tcsims[f"sor_po400"] = TimecourseSim(
            Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=200,
                changes={
                    **self.default_changes(),
                    "PODOSE_sor": Q_(400, "mg"),
                },
            )
        )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for dose in self.doses:
          mappings[f"fm_sor_po400"] = FitMapping(
            self,
            reference=FitData(
                self,
                dataset=f"SOF400",
                xid="time",
                yid="mean",
                yid_sd=None,
                count="count",
            ),
            observable=FitData(
                self, task=f"task_sor_po400", xid="time", yid="[Cve_sor]",
            ),
        )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        name = "Fig4"
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=1,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit="hour"),
            yaxis=Axis(self.label_sor, unit=self.units["[Cve_sor]"]),
            legend=True
        )

        # simulation
        for dose in self.doses:
          plots[0].add_data(
            task="task_sor_po400",
            xid="time",
            yid="[Cve_sor]",
            label="Sim sorafenib 400 [mg]",
            color=self.color_for_dose(dose),
        )

        # data
        for dose in self.doses:
          plots[0].add_data(
            dataset="SOF400",
            xid="time",
            yid="mean",
            yid_sd=None,
            count="count",
            label=f"400 [mg]",
            color=self.color_for_dose(dose),
        )

        return {
            fig.sid: fig
        }


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Andriamanana2013.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Andriamanana2013, output_dir="Andriamanana2013")
