from copy import deepcopy
from typing import Dict
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Axis, Figure
from pkdb_models.models import sorafenib
from sbmlsim.simulation import Timecourse, TimecourseSim
from pkdb_models.models.sorafenib.helpers import run_experiments
from pkdb_models.models.sorafenib.experiments.base_experiment import SorafenibSimulationExperiment


class Ferrario2016(SorafenibSimulationExperiment):
    """Simulation experiment of Ferrario2016."""
    doses = [200]
    substances = ["sor", "m2"]
    yids = ["[Cve_sor]", "[Cve_m2]"]
    groups = {
        "sor": ["SOF20_sor", "SOF22_sor", "SOF23_sor", "SOF27_sor", "SOF32_sor", "SOF33_sor"],
        "m2": ["SOF20_m2", "SOF22_m2", "SOF23_m2", "SOF27_m2", "SOF32_m2", "SOF33_m2"]
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig4"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)

            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)

                if label.endswith("_sor"):
                    dset.unit_conversion("value", 1 / self.Mr.sor)
                elif label.endswith("_m2"):
                    dset.unit_conversion("value", 1 / self.Mr.m2)
                dsets[label] = dset

        # print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        # multiple dosing
        Q_ = self.Q_
        tcsims = {}

        for dose in self.doses:
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

            for substance in self.substances:
                tcsims[f"{substance}_po_{dose}"] = TimecourseSim(
                    [tc0] + [deepcopy(tc1) for _ in range(13)],  # assumed 2 weeks of treatment,
                    time_offset=-13 * 24 * 60
                )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for dose in self.doses:
          for k, substance in enumerate(self.substances):
            datasets = self.groups[substance]
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
                        self, task=f"task_sor_po_{dose}", xid="time", yid=self.yids[k],
                    ),
                )

        return mappings


    @property
    def figures(self) -> Dict[str, Figure]:
        name = "Fig4"
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=2,
            num_cols=2,
            name=f"{self.__class__.__name__}",
        )
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="hour"), legend=False)
        plots[0].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])
        plots[1].set_yaxis(self.label_sor, unit=self.units["[Cve_sor]"])
        plots[2].set_yaxis(self.label_m2, unit=self.units["[Cve_m2]"])
        plots[3].set_yaxis(self.label_m2, unit=self.units["[Cve_m2]"])
        plots[1].legend = True
        plots[3].legend = True
        for k in [1, 3]:
            plots[k].xaxis.min = -12

        # Simulation
        for ks, substance in enumerate(self.substances):
            for k in [2*ks, 2*ks+1]:
                for dose in self.doses:
                   plots[k].add_data(
                    task=f"task_{substance}_po_{dose}",
                    xid="time",
                    yid=f"[Cve_{substance}]",
                    label=f"Sim {dose} mg {substance}",
                    color=self.color_for_dose(dose),
                )

        # Data

        for ks, substance in enumerate(self.substances):
            datasets = self.groups[substance]
            for dataset in datasets:
                for k in [2 * ks, 2 * ks + 1]:
                    plots[k].add_data(
                        dataset=dataset,
                        xid="time",
                        yid="value",
                        yid_sd=None,
                        count="count",
                        label=f"{dataset}",
                        color=self.color_for_dose(dose),
                    )
        return {
            fig.sid: fig
        }


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Ferrario2016.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Ferrario2016, output_dir="Ferrario2016")
