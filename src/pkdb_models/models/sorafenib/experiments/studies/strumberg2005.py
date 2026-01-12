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


class Strumberg2005(SorafenibSimulationExperiment):
    """Simulation experiment of Strumberg2005"""
    doses = [100, 200, 400, 600, 800]
    labels = ["SOF100_SD","SOF100_MD","SOF200_SD", "SOF200_MD","SOF400_SD", "SOF400_MD","SOF600_SD", "SOF600_MD", "SOF800_SD", "SOF800_MD"]


    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)

            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                dset.unit_conversion("mean", 1 / self.Mr.sor)
                dsets[f"sor_{label}"] = dset

        # print(dsets.keys())
        # print(dsets)

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for label in self.labels:
            for dose in self.doses:

                if label.endswith("SD"):
                    tcsims[f"sor_po_{dose}_SD"] = TimecourseSim(
                        Timecourse(
                            start=0,
                            end=12 * 60,  # [min]
                            steps=1000,
                            changes={
                                **self.default_changes(),
                                "PODOSE_sor": Q_(dose, "mg"),
                            },
                        )
                    )
                elif label.endswith("MD"):
                    tc0 = Timecourse(
                        start=0,
                        end=24 * 60,  # [min]
                        steps=200,
                        changes={
                            **self.default_changes(),
                            "PODOSE_sor": Q_(dose, "mg"),
                        },
                    )
                    tc1 = Timecourse(
                        start=0,
                        end=24 * 60,  # [min]
                        steps=200,
                        changes={
                            "PODOSE_sor": Q_(dose, "mg"),
                        },
                    )

                    tcsims[f"sor_po_{dose}_MD"] = TimecourseSim(
                        [tc0] + [deepcopy(tc1) for _ in range(14)],  # assumed 2 weeks of treatment
                         time_offset=-14 * 24 * 60
                    )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        datainfo = {
            100: ["SOF100_SD", "SOF100_MD"],
            200: ["SOF200_SD", "SOF200_MD"],
            400: ["SOF400_SD", "SOF400_MD"],
            600: ["SOF600_SD", "SOF600_MD"],
            800: ["SOF800_SD", "SOF800_MD"],
        }

        for dose in self.doses:
            datasets = datainfo[dose]
            for dataset in datasets:
                dosing = dataset.split("_")[-1]
                mappings[f"fm_sor_po_{dose}_{dosing}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"sor_{dataset}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_sor_po_{dose}_{dosing}", xid="time", yid="[Cve_sor]",
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:

        name = "Fig1"  #Fig2
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=2,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="hour"), legend=False)
        plots[0].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])
        plots[1].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])

        # simulation
        for dose in self.doses:
            for k, dosing in enumerate(["SD", "MD"]):
                plots[k].add_data(
                    task=f"task_sor_po_{dose}_{dosing}",
                    xid="time",
                    yid="[Cve_sor]",
                    label=f"Sim {dose} [mg])",
                    color=self.color_for_dose(dose),
                )

        # Data
        for dose in self.doses:
            plots[0].add_data(
                dataset=f"sor_SOF{dose}_SD",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose}[mg]",
                color=self.color_for_dose(dose),
            )
            plots[1].add_data(
                dataset=f"sor_SOF{dose}_MD",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{dose}[mg]",
                color=self.color_for_dose(dose),
            )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Strumberg2005.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Strumberg2005, output_dir="Strumberg2005")
