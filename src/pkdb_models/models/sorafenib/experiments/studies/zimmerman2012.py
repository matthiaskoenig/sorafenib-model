from copy import deepcopy
from typing import Dict
from pkdb_models.models import sorafenib
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from pkdb_models.models.sorafenib.experiments.base_experiment import SorafenibSimulationExperiment
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.sorafenib.helpers import run_experiments


class Zimmerman2012(SorafenibSimulationExperiment):
    """Simulation experiment of Zimmerman2012."""
    doses = [150, 200]
    substances = ["sor", "m2", "sg"]
    yids = ["[Cve_sor]", "[Cve_sg]", "[Cve_sg]"]
    datainfo = {
        "sor": ["SOF150_sor", "SOF200_sor"],
        "m2": ["SOF150_m2", "SOF200_m2"],
        "sg": ["SOF150_sg", "SOF200_sg"],
                }

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
                elif label.endswith("_sg"):
                    dset.unit_conversion("mean", 1 / self.Mr.sg)
                dsets[label] = dset

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

            tcsims[f"sor_po{dose}"] = TimecourseSim(
                [tc_first_day] + [deepcopy(tc_dosing) for _ in range(4)],
                time_offset=-4 * 24 * 60
            )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for dose in self.doses:
            for k, substance in enumerate(self.substances):
                datasets = self.datainfo[substance]
                for dataset in datasets:
                    mappings[f"fm_sor_po_{dataset}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=dataset,
                            xid="time",
                            yid="mean",
                            yid_sd=None,
                            count="count",
                        ),
                        observable=FitData(
                            self, task=f"task_sor_po{dose}", xid="time", yid=self.yids[k],
                        ),
                    )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        name = "Fig1"
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=3,
            num_cols=1,
            name=f"{self.__class__.__name__}",
        )

        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="hour"), legend=True)
        plots[0].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])
        plots[1].set_yaxis(self.label_m2, unit=self.units["[Cve_m2]"])
        plots[2].set_yaxis(self.label_sg, unit=self.units["[Cve_sg]"])

        # Simulation
        for k, substance in enumerate(self.substances):
            datasets = self.datainfo[substance]
            for dose in self.doses:
                plots[k].add_data(
                    task=f"task_sor_po{dose}",
                    xid="time",
                    yid=self.yids[k],
                    label=f"Sim {dose} mg {substance}",
                    color=self.color_for_dose(dose),
                )

            for j, dataset in enumerate(datasets):
                plots[k].add_data(
                    dataset=dataset,
                    xid="time",
                    yid="mean",
                    count="count",
                    yid_sd=None,
                    label=f"{substance} [mg] ({dataset})",
                    color=self.color_for_dose(dose),
                 )

        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Zimmerman2012.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Zimmerman2012, output_dir="Zimmerman2012")
