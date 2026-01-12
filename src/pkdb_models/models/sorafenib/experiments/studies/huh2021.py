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


class Huh2021(SorafenibSimulationExperiment):
    """Simulation experiment of Huh2021.

    # 10 µM = 10 µmol/l => 10 µmol/l * Mr = 10 µmol/l * 464.826 g/mole = 4648 µg/l
    """
    doses = [100, 150, 200]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
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
        for dose in self.doses:
            tcsims[f"sor_po_{dose}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=170 * 60,  # [min]
                    steps=1000,
                    changes={
                        **self.default_changes(),
                        "PODOSE_sor": Q_(dose, "mg"),
                    },
                )
            )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        datainfo = {
            100: ["SYO100_P1", "SYO100_P2"],
            150: ["SYO150_P1", "SYO150_P2"],
            200: ["SYO200_P1", "SOF200_P1", "SYO200_P2", "SOF200_P2"]
        }
        for dose in self.doses:
            datasets = datainfo[dose]
            for dataset in datasets:
                mappings[f"fm_sor_po_{dose}_{dataset}"] = FitMapping(
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
                        self, task=f"task_sor_po_{dose}", xid="time", yid="[Cve_sor]",
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
            plots[0].add_data(
                task=f"task_sor_po_{dose}",
                xid="time",
                yid="[Cve_sor]",
                label=f"Sim {dose} [mg]",
                color=self.color_for_dose(dose),
            )

        # data
        datainfo = {
            100: ["SYO100_P1", "SYO100_P2"],
            150: ["SYO150_P1", "SYO150_P2"],
            200: ["SYO200_P1", "SOF200_P1", "SYO200_P2", "SOF200_P2"]
        }

        for dose in self.doses:
            datasets = datainfo[dose]
            for dataset in datasets:
                plots[0].add_data(
                    dataset=f"sor_{dataset}",
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
    out = sorafenib.RESULTS_PATH_SIMULATION / Huh2021.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Huh2021, output_dir="Huh2021")
