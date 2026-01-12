
from copy import deepcopy
from typing import Dict

import numpy as np
from matplotlib import pyplot as plt
from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.plot import Axis, Figure
from pkdb_models.models import sorafenib
from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.simulation import Timecourse, TimecourseSim
from sbmlutils.console import console

from pkdb_models.models.sorafenib.experiments.base_experiment import SorafenibSimulationExperiment
from pkdb_models.models.sorafenib.helpers import run_experiments


class Hornecker2012(SorafenibSimulationExperiment):
    """Simulation experiment of Hornecker2012.

    TODO: calculate pharmacokinetic parameters and plot tmax ~ dose, auc_end ~ dose;
    also bioavailability data?
    """
    doses = [400, 800, 1200, 1600, 2000, 2400]
    n_doses = 15

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for mtype, df_label in df.groupby("measurement_type"):
                dset = DataSet.from_df(df_label, self.ureg)

                dsets[f"{mtype}"] = dset
        # print(dsets.keys())
        # print(dsets)
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

            tcsims[f"sor_po{dose}"] = TimecourseSim(
                [tc0] + [deepcopy(tc1) for _ in range(self.n_doses-1)],  # assumed 2 weeks of treatment,
                # time_offset=-14 * 24 * 60
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:
        name = "Fig3"
        fig = Figure(
            experiment=self,
            sid=f"{name}_PO_plasma",
            num_rows=1,
            num_cols=2,
            name=None,
        )

        plots = fig.create_plots(xaxis=Axis(self.label_time, unit="day"), legend=True)
        plots[0].set_yaxis(self.labels["[Cve_sor]"], unit=self.units["[Cve_sor]"])
        plots[1].set_yaxis(self.labels["[Cve_sor]"], unit=self.units["[Cve_sor]"])
        plots[0].xaxis.min = 0
        plots[0].xaxis.max = 1

        #Simulation
        for k, dose in enumerate(self.doses):
            for kplot in [0, 1]:
                plots[kplot].add_data(
                    task=f"task_sor_po{dose}",
                    xid="time",
                    yid="[Cve_sor]",
                    label=f"Sim {dose} mg",
                    color=self.color_for_dose(dose),
                )

        return {
            fig.sid: fig
        }

    def figures_mpl(self) -> Dict[str, FigureMPL]:

        # calculate the pharmacokinetic parameters for the last dose
        self.pk_dfs_multi = self.calculate_sorafenib_pk(
            tstart=(self.n_doses-1) * 24 * 60,  # shift time slot to last day
            tend=self.n_doses * 24 * 60,  # shift time slot to last day
        )
        self.pk_dfs_single = self.calculate_sorafenib_pk(
            tstart=0.0,  # shift time slot to last day
            tend=24 * 60,  # shift time slot to last day
        )

        return {
           **self.figures_mpl_pharmacokinetics(),
        }

    def figures_mpl_pharmacokinetics(self):
        """Visualize dependency of pharmacokinetics parameters."""
        #global df
        Q_ = self.Q_
        parameters = [
            "auc",
            "tmax",
            "cmax",
            "kel",
        ]

        f, axes = plt.subplots(nrows=1, ncols=4, figsize=(24, 6))
        f.subplots_adjust(wspace=0.3)
        axes = axes.flatten()

        doses = np.array(self.doses)
        for dosing in ["single", "multi"]:
            if dosing == "single":
                pk_dfs = self.pk_dfs_single
            elif dosing == "multi":
                pk_dfs = self.pk_dfs_multi

            for k_key, pk_key in enumerate(parameters):
                ax = axes[k_key]

                # collect the PK parameter vector
                pk_vec = np.zeros_like(doses, dtype=float)
                for k_dose, dose in enumerate(doses):
                    df = pk_dfs[f"sor_po{dose}"]
                    pk_vec[k_dose] = df[pk_key].values[0]

                x = Q_(doses, "mg")
                y = Q_(pk_vec, df[f"{pk_key}_unit"][0])
                y = y.to(self.pk_units[pk_key])

                ax.plot(
                    x,
                    y,
                    marker="o" if dosing == "single" else "^",
                    linestyle="-",
                    markeredgecolor="black",
                    markersize=8,
                    label=f"simulation ({dosing})",
                )

                ax.set_xlabel("Dose [mg]", fontdict=SorafenibSimulationExperiment.font)
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )
                ax.set_ylim(bottom=0.0)
                # ax.legend(fontsize=SorafenibSimulationExperiment.legend_font_size)


        # Plot experimental data
        df_auc = self._datasets["auc_end"]
        doses = df_auc["dose"].values
        axes[0].errorbar(
            x=doses,
            y=df_auc["median"]/self.Mr.sor.magnitude,  # [mg/l/hr] -> [mmole/l/hr]
            yerr=[
                np.abs(df_auc["median"]-df_auc["wlow"])/self.Mr.sor.magnitude,
                np.abs(df_auc["median"]-df_auc["wup"])/self.Mr.sor.magnitude,
            ],
            color="black",
            marker="s",
            linestyle="-",
            markeredgecolor="black",
            label="Hornecker2012 (median)",
            markersize=10,
        )

        df_tmax = self._datasets["tmax"]
        axes[1].errorbar(
            x=doses,
            y=df_tmax["median"],  # [mg/l/hr] -> [mmole/l/hr]
            yerr=[
                np.abs(df_tmax["median"]-df_tmax["wlow"]),
                np.abs(df_tmax["median"]-df_tmax["wup"]),
            ],
            color="black",
            marker="s",
            linestyle="-",
            markeredgecolor="black",
            label="Hornecker2012 (median)",
            markersize=10,
        )

        for ax in axes:
            ax.legend()
            ax.set_xlim(left=0.0)

        axes[0].set_ylim(top=1.0)
        axes[1].set_ylim(top=10.0)
        axes[2].set_ylim(top=50)
        axes[3].set_ylim(top=0.05)

        return {
            "fig_pk": f
        }


if __name__ == "__main__":
    out = sorafenib.RESULTS_PATH_SIMULATION / Hornecker2012.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Hornecker2012, output_dir="Hornecker2012")
