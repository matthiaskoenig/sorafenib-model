from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console
from pkdb_models.models import sorafenib

from pkdb_models.models.sorafenib.experiments.base_experiment import SorafenibSimulationExperiment
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.sorafenib.helpers import run_experiments


class Mammatas2020(SorafenibSimulationExperiment):
    """Simulation experiment of Mammatas2020."""
    doses = [1000, 2000, 2400, 2800]
    groups = ["SOF_1000", "SOF_2000", "SOF_2400", "SOF_2800"]

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)

            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                # console.print(dset.to_string())
                dset.unit_conversion("mean", 1 / self.Mr.sor)
                dsets[f"sor_{label}"] = dset

        # print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for dose in self.doses:
            tcsims[f"sor_po_{dose}"] = TimecourseSim(
                Timecourse(
                    start=0,
                    end=72 * 60,
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
        for dose in self.doses:
            for group in self.groups:
                mappings[f"fm_sor_po_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"sor_{group}",
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
        name = "Fig1"
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
            dataset = f"sor_SOF_{dose}"
            plots[0].add_data(
                task=f"task_sor_po_{dose}",
                xid="time",
                yid="[Cve_sor]",
                label=f"Sim {dose} [mg]",
                color=self.color_for_dose(dose),
            )

            # Adding data
            for dose in self.doses:
                plots[0].add_data(
                    dataset=dataset,
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=f"{dose} [mg]",
                    color=self.color_for_dose(dose),
                )
        return {fig.sid: fig}


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Mammatas2020.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Mammatas2020, output_dir="Mammatas2020")
