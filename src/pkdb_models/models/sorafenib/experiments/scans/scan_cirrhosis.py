from typing import Dict
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.cm as cm

from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension

from sbmlsim.plot.serialization_matplotlib import FigureMPL
from sbmlsim.plot.serialization_matplotlib import plt
from pkdb_models.models.sorafenib.experiments.base_experiment import SorafenibSimulationExperiment
from pkdb_models.models.sorafenib.helpers import run_experiments


class CirrhosisScan(SorafenibSimulationExperiment):
    """Scan the effect of hepatic function on sorafenib pharmacokinetics."""

    scan_map = {
        "hepatic": {
            "parameter": "f_cirrhosis",
            "range": np.linspace(0, 0.9, num=19),
            "units": "dimensionless",
            "label": "Cirrhosis [-]"
        },
    }

    def simulations(self) -> Dict[str, ScanSim]:
        Q_ = self.Q_
        tcscans = {}
        for scan_key, scan_data in self.scan_map.items():
            for renal_key in self.renal_map.keys():

                # FIXME: Only calculation on single dose! This has to be performed on multiple doses
                tcscans[f"scan_{scan_key}_{renal_key}"] = ScanSim(
                    simulation=TimecourseSim(
                        Timecourse(
                            start=0,
                            end=24 * 60,
                            steps=200,
                            changes={
                                **self.default_changes(),
                                "PODOSE_sor": Q_(400, 'mg'),
                                "KI__f_renal_function": Q_(self.renal_map[renal_key], "dimensionless"),
                            },
                        )
                    ),
                    dimensions=[
                        Dimension(
                            "dim_scan",
                            changes={
                                scan_data["parameter"]: Q_(
                                    scan_data["range"], scan_data["units"]
                                )
                            },
                        ),
                    ],
                )

        return tcscans

    def figures_mpl(self) -> Dict[str, FigureMPL]:
        # calculate the pharmacokinetic parameters
        self.pk_dfs = self.calculate_sorafenib_pk(

        )

        return {
           **self.figures_mpl_timecourses(),
           **self.figures_mpl_pharmacokinetics(),
        }

    def figures_mpl_timecourses(self) -> Dict[str, FigureMPL]:
        """Timecourse plots for key variables depending on renal impairment degree."""
        sids = [
            "[Cve_sor]",
            # "Aurine_pra",
            # "[LI__sor]",
            # "LI__pra_bi",
            # "Afeces_pra",
            "Afeces_sor",
            "Aurine_sg",
            "Afeces_sg",
        ]
        n_sids = len(sids)
        cmap = matplotlib.cm.get_cmap('seismic')
        figures = {}
        for scan_key, scan_data in self.scan_map.items():
            f, axes = plt.subplots(nrows=n_sids, ncols=len(self.renal_map),
                                   figsize=(5 * len(self.renal_map), 5 * len(sids)),
                                   sharey="row", sharex="col"
                                )
            f.subplots_adjust(wspace=0.1, hspace=0.1)
            f.suptitle(f"Scan {scan_data['label']}", fontsize=self.suptitle_font_size)

            ymax = {}
            for krow, sid in enumerate(sids):
                ymax[sid] = 0.0
                for kcol, renal_key in enumerate(self.renal_map.keys()):
                    ax = axes[krow, kcol]
                    # get data
                    Q_ = self.Q_
                    xres = self.results[f"task_scan_{scan_key}_{renal_key}"]

                    # scanned dimension
                    scandim = xres._redop_dims()[0]

                    parameter_id = scan_data["parameter"]
                    par_vec = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])
                    t_vec = xres.dim_mean("time").to(self.units["time"])

                    # update ymax
                    y = Q_(
                        xres[sid].values, xres.uinfo[sid]
                    ).to(self.units[sid])
                    if np.max(y.magnitude) > ymax[sid]:
                        ymax[sid] = np.max(y.magnitude)

                    for k_par, par in enumerate(par_vec):
                        c_vec = Q_(
                            xres[sid].sel({scandim: k_par}).values, xres.uinfo[sid]
                        ).to(self.units[sid])

                        # plot all curves for the scan
                        facecolor = self.renal_colors[renal_key]
                        if facecolor != "black":
                            for axis in ['top', 'bottom', 'left', 'right']:
                                ax.spines[axis].set_linewidth(4)
                                ax.spines[axis].set_color(facecolor)

                        # 0.0 - 0.9
                        linewidth = 2.0
                        if np.isclose(0.0, par):
                            color = "black"
                        elif par > 0.0:
                            cvalue = 0.5 - par/0.9 * 0.5
                            color = cmap(cvalue)

                        ax.plot(t_vec.magnitude, c_vec.magnitude, color=color, linewidth=linewidth)

                    ax: matplotlib.axes.Axes
                    if krow == 0:
                        ax.set_title(renal_key, fontdict={
                            'fontsize': SorafenibSimulationExperiment.tick_font_size,
                            'fontweight': 'bold'
                        })

                    if krow == 4:
                        ax.set_xlabel(f"{self.label_time} [{self.units['time']}]",
                                      fontdict=SorafenibSimulationExperiment.font)

                    ax.tick_params(axis="x", labelsize=SorafenibSimulationExperiment.tick_font_size)
                    ax.tick_params(axis="y", labelsize=SorafenibSimulationExperiment.tick_font_size)



                    if kcol == 0:
                        ax.set_ylabel(f"{self.labels[sid]} [{self.units[sid]}]",
                                      fontdict=SorafenibSimulationExperiment.font)


            for krow, sid in enumerate(sids):
                for kcol, renal_key in enumerate(self.renal_map.keys()):
                    ax = axes[krow, kcol]
                    ax.set_ylim(bottom=0, top=1.2 * ymax[sid])

            for kcol in range(4):
                ax = axes[0, kcol]
                ax.set_ylim(bottom=0, top=4)

                # add colorbar
                cb_ax = f.add_axes([0.92, 0.2, 0.02, 0.6])
                norm = matplotlib.colors.Normalize(vmin=0, vmax=2, clip=False)

                cbar = f.colorbar(cm.ScalarMappable(norm=norm, cmap='seismic_r'), cax=cb_ax)
                cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

                cbar.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0],
                                    **{"size": 30, "weight": "medium"}
                                    )
                cbar.ax.set_ylabel(scan_data["label"], rotation=270, **{"size": 25, "weight": "bold"})

        figures[f"fig_timecourse_{scan_key}"] = f
        return figures


    def figures_mpl_pharmacokinetics(self):
        """Visualize dependency of pharmacokinetics parameters."""
        Q_ = self.Q_
        figures = {}
        parameters = [
            "aucinf",
            "kel",
            "tmax",
            "thalf",
            "cl",
            # "cl_hepatic",
            # "cl_renal",
        ]

        for scan_key, scan_data in self.scan_map.items():
            f, axes = plt.subplots(nrows=1, ncols=5, figsize=(6*5, 6*1)) #6*3, 6*2
            f.subplots_adjust(wspace=0.3)
            axes = axes.flatten()
            for k, pk_key in enumerate(parameters):
                ax = axes[k]

                # vertical lines
                cirrhosis_labels = ["Healthy", "CTP A", "CPT B", "CPT C"]
                cirrhosis_colors = list(self.cirrhosis_colors.values())
                for k, value in enumerate(self.cirrhosis_map.values()):
                    ax.axvline(x=value, color=cirrhosis_colors[k], linestyle="--")
                    xy = (value-0.04, 1.02)
                    if k == 3:
                        xy = (value, 1.02)
                    ax.annotate(
                        text=cirrhosis_labels[k],
                        xy=xy,
                        xycoords=("data", "axes fraction"),
                        fontweight="bold",
                        # color=cirrhosis_colors[k]
                    )

                ax.tick_params(axis="x", labelsize=SorafenibSimulationExperiment.tick_font_size)
                ax.tick_params(axis="y", labelsize=SorafenibSimulationExperiment.tick_font_size)

                for renal_key in self.renal_map:
                    sim_key = f"scan_{scan_key}_{renal_key}"
                    xres = self.results[f"task_{sim_key}"]
                    df = self.pk_dfs[sim_key]

                    # This was scanned
                    parameter_id = scan_data["parameter"]
                    x_vec = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])

                    pk_vec = df[f"{pk_key}"]
                    pk_vec = pk_vec.to_numpy()

                    x = x_vec
                    y = Q_(pk_vec, df[f"{pk_key}_unit"][0])
                    y = y.to(self.pk_units[pk_key])

                    ax.plot(
                        x,
                        y,
                        marker="o",
                        linestyle="-",
                        color=self.renal_colors[renal_key],
                        markeredgecolor="black",
                        label=renal_key,
                    )

                # ax.tick_params(axis="x", labelsize=SorafenibSimulationExperiment.tick_font_size)
                # ax.tick_params(axis="y", labelsize=SorafenibSimulationExperiment.tick_font_size)
                ax.set_xlabel(scan_data["label"], fontdict=SorafenibSimulationExperiment.font)
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )
                ax.set_ylim(bottom=0.0)
                ax.legend(fontsize=SorafenibSimulationExperiment.legend_font_size)
                ax.legend()

            figures[f"fig_{scan_key}"] = f
        return figures


if __name__ == "__main__":
    run_experiments(CirrhosisScan, output_dir="scan_hepatic")
