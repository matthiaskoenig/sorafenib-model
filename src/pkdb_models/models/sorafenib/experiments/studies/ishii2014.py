from copy import deepcopy
from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models.sorafenib.experiments.base_experiment import SorafenibSimulationExperiment
from sbmlsim.plot import Axis, Figure
from pkdb_models.models import sorafenib
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.sorafenib.helpers import run_experiments


class Ishii2014(SorafenibSimulationExperiment):
    """Simulation experiment of Ishii2014.

    FIXME: Hemodialysis !!!
    """
    doses =[200]
    substances = ["sor", "m2"]
    yids = ["[Cve_sor]", "[Cve_m2]"]
    labels = ["SOF_sor", "SOF_m2"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.endswith("_sor"):
                    dset.unit_conversion("mean", 1 / self.Mr.sor)
                elif label.endswith("_m2"):
                    dset.unit_conversion("mean", 1 / self.Mr.m2)
                dsets[label] = dset

        # print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for dose in self.doses:
          tc0 = Timecourse(
            start=0,
            end=24 * 60,  # [min]
            steps=200,
            changes={
                **self.default_changes(),
                "PODOSE_sor": Q_(200, "mg"),
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

        tcsims[f"sor_po_200"] = TimecourseSim(
            [tc0] + [deepcopy(tc1) for _ in range(8)],  # assumed 2 weeks of treatment,
            time_offset=-7 * 24 * 60
        )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for dose in self.doses:
          for k, substance in enumerate(self.substances):
            for label in self.labels:
                mappings[f"fm_sor_po_{label}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=label,
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_sor_po_200", xid="time", yid=self.yids[k],
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        name = "Fig1"
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

        #colors = {
            #"sor": "tab:green",
            #"m2": "tab:red"
        #}
        # Simulation
        for k, substance in enumerate(self.substances):
          for dose in self.doses:
            plots[k].add_data(
                task=f"task_sor_po_200",
                xid="time",
                yid=f"[Cve_{substance}]",
                label=f"Sim 200 [mg] ({substance})",
                color=self.color_for_dose(dose),
            )

            # data
            dataset = f"SOF_{substance}"
            for dose in self.doses:
              plots[k].add_data(
                dataset=dataset,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f" 200 [mg] ({dataset})",
                color=self.color_for_dose(dose),
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Ishii2014.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Ishii2014, output_dir="Ishii2014")
