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

class Huang2017(SorafenibSimulationExperiment):
    """Simulation experiment of Huang2017."""
    doses = [200, 400]
    substances = ["sor", "m2"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3", "Fig4"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)

            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                if label.endswith("_sor"):
                    dset.unit_conversion("mean", 1 / self.Mr.sor)
                elif label.endswith("_m2"):
                    dset.unit_conversion("mean", 1 / self.Mr.m2)
                dsets[label] = dset

        # print(dsets.keys())
        # print(dsets)

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for dose in self.doses:
          tcsims[f"sor_po_400"] = TimecourseSim(
            Timecourse(
                start=0,
                end=96 * 60,  # [min]
                steps=1000,
                changes={
                    **self.default_changes(),
                    "PODOSE_sor": Q_(400, "mg"),
                },
            )
        )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for substance in ["sor", "m2"]:
            mappings[f"fm_sor_po_{substance}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"SOF400_{substance}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_sor_po_400", xid="time", yid=f"[Cve_{substance}]"
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:

        name = "Fig3_4"
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=2,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="hour"), legend=True)
        plots[0].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])
        plots[1].set_yaxis(self.label_m2, unit=self.units["[Cve_m2]"])

        # simulation
        for k, substance in enumerate(self.substances):
            # simulation
          for dose in self.doses:
            plots[k].add_data(
                task=f"task_sor_po_400",
                xid="time",
                yid=f"[Cve_{substance}]",
                label=f"Sim {substance} 400 [mg])",
                color=self.color_for_dose(dose),
            )

            # data
            dataset = f"SOF400_{substance}"
            for dose in self.doses:
              plots[k].add_data(
                dataset=dataset,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"400 [mg] ({dataset})",
                color=self.color_for_dose(dose),
            )

        return {
            fig.sid: fig
        }


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Huang2017.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Huang2017, output_dir="Hornecker2012")
