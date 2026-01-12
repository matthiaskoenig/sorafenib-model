from typing import Dict

import matplotlib.axes
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import seaborn as sns
from sbmlsim.data import DataSet, load_pkdb_dataframe

from sbmlsim.simulation import Timecourse, TimecourseSim, ScanSim, Dimension
from sbmlsim.plot import Figure, Axis
from sbmlsim.plot.serialization_matplotlib import FigureMPL, MatplotlibFigureSerializer
from sbmlsim.plot.serialization_matplotlib import plt
from pkdb_models.models.sorafenib.experiments.base_experiment import SorafenibSimulationExperiment
from pkdb_models.models.sorafenib.helpers import run_experiments


class ParametersScan(SorafenibSimulationExperiment):
    """Scan the effect of renal function on sorafenib pharmacokinetics."""

    # FIXME: BETTER SCANS: "range": np.sort(np.append(np.logspace(-1, 1, num=num_points), [1.0])),  # [10^-1=0.1, 10^1=10]
    num_points = 19
    scan_map = {
        "renal": {
            "parameter": "KI__f_renal_function",
            "range": np.logspace(-1, 1, num=num_points),
            "units": "dimensionless",
            "label": "Renal function [-]"
        },
        # "mrp2": {
        #     "parameter": "LI__f_MRP2",
        #     "range": np.logspace(-1, 1, num=num_points),
        #     "units": "dimensionless",
        #     "label": "MRP2 (ABCC2) activity [-]"
        # },
    }

    def simulations(self) -> Dict[str, ScanSim]:
        Q_ = self.Q_
        tcscans = {}

        for scan_key, scan_data in self.scan_map.items():
            print(scan_data["range"])
            for cirrhosis_key in self.cirrhosis_map.keys():
                for renal_key in self.renal_map.keys():
                    tcscans[f"scan_{scan_key}_{cirrhosis_key}_{renal_key}"] = ScanSim(
                        simulation=TimecourseSim(
                            Timecourse(
                                start=0,
                                end=24 * 60,
                                steps=200,
                                changes={
                                    **self.default_changes(),
                                    "PODOSE_sor": Q_(400, 'mg'),
                                    "f_cirrhosis": Q_(self.cirrhosis_map[cirrhosis_key], "dimensionless"),
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
        self.pk_dfs = self.calculate_sorafenib_pk()

        return {
            **self.figures_mpl_timecourses_cirrhosis(),
            **self.figures_mpl_timecourses_renal(),
            **self.figures_mpl_pharmacokinetics(),
        }

    def figures_mpl_timecourses_cirrhosis(self) -> Dict[str, FigureMPL]:
        """Timecourse plots for key variables depending on cirrhosis degree."""
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
        cmap = matplotlib.cm.get_cmap('seismic')
        figures = {}
        for scan_key, scan_data in self.scan_map.items():
            f, axes = plt.subplots(nrows=len(sids), ncols=len(self.cirrhosis_map), figsize=(5 * len(self.cirrhosis_map), 5 * len(sids)),
                                   sharex="col", sharey="row")
            f.subplots_adjust(wspace=0.05, hspace=0.05)
            f.suptitle(f"Scan {scan_data['label']} in cirrhosis", fontsize=self.suptitle_font_size)

            ymax = {}
            for krow, sid in enumerate(sids):
                ymax[sid] = 0.0
                for kcol, cirrhosis_key in enumerate(self.cirrhosis_map.keys()):

                    ax = axes[krow, kcol]
                    # get data
                    Q_ = self.Q_
                    xres = self.results[f"task_scan_{scan_key}_{cirrhosis_key}_{list(self.renal_map.keys())[0]}"]

                    # scanned dimension
                    scandim = xres._redop_dims()[0]

                    parameter_id = scan_data["parameter"]
                    par_vec = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])
                    print(par_vec)
                    t_vec = xres.dim_mean("time").to(self.units["time"])

                    for k_par, par in enumerate(par_vec):
                        c_vec = Q_(
                            xres[sid].sel({scandim: k_par}).values, xres.uinfo[sid]
                        ).to(self.units[sid])

                        # update ymax
                        if np.max(c_vec.magnitude) > ymax[sid]:
                            ymax[sid] = np.max(c_vec.magnitude)

                        # plot all curves for the scan
                        facecolor = self.cirrhosis_colors[cirrhosis_key]
                        if facecolor != "black":
                            for axis in ['top', 'bottom', 'left', 'right']:
                                ax.spines[axis].set_linewidth(4)
                                ax.spines[axis].set_color(facecolor)

                        # 0.1 - 1.9
                        linewidth = 2.0
                        if np.isclose(1.0, par):
                            color = "black"
                            t_vec_one = t_vec
                            c_vec_one = c_vec
                        else:
                            cvalue = 1 - ((par - 0.1) / 1.8)  # red less function, blue more function
                            color = cmap(cvalue)

                        ax.plot(t_vec.magnitude, c_vec.magnitude, color=color, linewidth=linewidth)

                    ax.plot(t_vec_one.magnitude, c_vec_one.magnitude, color="black", linewidth=2.0)

                    if krow == 0:
                        ax.set_title(cirrhosis_key, fontdict={
                            'fontsize': SorafenibSimulationExperiment.tick_font_size,
                            'fontweight': 'bold'
                        })

                    ax.tick_params(axis="x", labelsize=SorafenibSimulationExperiment.tick_font_size)
                    ax.tick_params(axis="y", labelsize=SorafenibSimulationExperiment.tick_font_size)

                    if krow == 4:
                        ax.set_xlabel(f"{self.label_time} [{self.units['time']}]", fontdict=SorafenibSimulationExperiment.font)
                    else:
                        ax.set_xticks(np.linspace(0, 30, 6))

                    if kcol == 0:
                        ax.set_ylabel(f"{self.labels[sid]} [{self.units[sid]}]",
                                  fontdict=SorafenibSimulationExperiment.font)
                    ax.set_ylim(bottom=0.0)

            for krow, sid in enumerate(sids):
                for kcol, cirrhosis_key in enumerate(self.cirrhosis_map.keys()):
                    ax = axes[krow, kcol]
                    ax.set_ylim(bottom=0, top=1.05 * ymax[sid])

            # add colorbar
            cb_ax = f.add_axes([0.92, 0.2, 0.02, 0.6])
            norm = matplotlib.colors.Normalize(vmin=0, vmax=2, clip=False)

            cbar = f.colorbar(cm.ScalarMappable(norm=norm, cmap='seismic_r'), cax=cb_ax)
            cbar.set_ticks([0, 0.5, 1, 1.5, 2])

            cbar.set_ticklabels([0, 0.5, 1, 1.5, 2],
                                 **{"size": 30, "weight": "medium"}
                                 )
            cbar.ax.set_ylabel(scan_data["label"], rotation=270, **{"size": 25, "weight": "bold"})


            figures[f"fig_timecourse_cirrhosis_{scan_key}"] = f
        return figures

    def figures_mpl_timecourses_renal(self) -> Dict[str, FigureMPL]:
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
        cmap = matplotlib.cm.get_cmap('seismic')
        figures = {}
        for scan_key, scan_data in self.scan_map.items():
            f, axes = plt.subplots(nrows=len(sids), ncols=len(self.renal_map),
                                   figsize=(6 * len(self.renal_map), 6 * len(sids)),
                                   sharex="col", sharey="row")
            f.subplots_adjust(wspace=0.02, hspace=0.02)
            f.suptitle(f"Scan {scan_data['label']} in renal impairment", fontsize=self.suptitle_font_size)

            ymax = {}
            for krow, sid in enumerate(sids):
                ymax[sid] = 0.0
                for kcol, renal_key in enumerate(self.renal_map.keys()):

                    ax = axes[krow, kcol]
                    # get data
                    Q_ = self.Q_
                    xres = self.results[f"task_scan_{scan_key}_{list(self.cirrhosis_map.keys())[0]}_{renal_key}"]

                    # scanned dimension
                    scandim = xres._redop_dims()[0]

                    parameter_id = scan_data["parameter"]
                    par_vec = Q_(xres[parameter_id].values[0], xres.uinfo[parameter_id])
                    t_vec = xres.dim_mean("time").to(self.units["time"])

                    for k_par, par in enumerate(par_vec):
                        c_vec = Q_(
                            xres[sid].sel({scandim: k_par}).values, xres.uinfo[sid]
                        ).to(self.units[sid])

                        # update ymax
                        if np.max(c_vec.magnitude) > ymax[sid]:
                            ymax[sid] = np.max(c_vec.magnitude)

                        # plot all curves for the scan
                        facecolor = self.renal_colors[renal_key]
                        if facecolor != "black":
                            for axis in ['top', 'bottom', 'left', 'right']:
                                ax.spines[axis].set_linewidth(4)
                                ax.spines[axis].set_color(facecolor)

                        # 0.1 - 1.9
                        linewidth = 2.0
                        if np.isclose(1.0, par):
                            color = "black"
                            t_vec_one = t_vec
                            c_vec_one = c_vec
                        else:
                            cvalue = 1 - ((par - 0.1) / 1.8)  # red less function, blue more function
                            color = cmap(cvalue)

                        ax.plot(t_vec.magnitude, c_vec.magnitude, color=color, linewidth=linewidth)

                    ax.plot(t_vec_one.magnitude, c_vec_one.magnitude, color="black", linewidth=2.0)

                    if krow == 0:
                        ax.set_title(renal_key, fontdict={
                            'fontsize': SorafenibSimulationExperiment.tick_font_size,
                            'fontweight': 'bold'
                        })

                    ax.tick_params(axis="x", labelsize=SorafenibSimulationExperiment.tick_font_size)
                    ax.tick_params(axis="y", labelsize=SorafenibSimulationExperiment.tick_font_size)

                    if krow == 4:
                        ax.set_xlabel(f"{self.label_time} [{self.units['time']}]", fontdict=SorafenibSimulationExperiment.font)
                    else:
                        ax.set_xticks(np.linspace(0, 30, 6))

                    if kcol == 0:
                        ax.set_ylabel(f"{self.labels[sid]} [{self.units[sid]}]",
                                  fontdict=SorafenibSimulationExperiment.font)
                    ax.set_ylim(bottom=0.0)

            for krow, sid in enumerate(sids):
                for kcol, renal_key in enumerate(self.renal_map.keys()):
                    ax = axes[krow, kcol]
                    ax.set_ylim(top=1.05 * ymax[sid])

            # add colorbar
            cb_ax = f.add_axes([0.92, 0.2, 0.02, 0.6])
            norm = matplotlib.colors.Normalize(vmin=0, vmax=2, clip=False)

            cbar = f.colorbar(cm.ScalarMappable(norm=norm, cmap='seismic_r'), cax=cb_ax)
            cbar.set_ticks([0, 0.5, 1, 1.5, 2])

            cbar.set_ticklabels([0, 0.5, 1, 1.5, 2],
                                **{"size": 30, "weight": "medium"}
                                )
            cbar.ax.set_ylabel(scan_data["label"], rotation=270, **{"size": 25, "weight": "bold"})

            figures[f"fig_timecourse_renal_{scan_key}"] = f
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
            f, axes = plt.subplots(nrows=1, ncols=5, figsize=(6 * 5, 6 * 1))
            f.subplots_adjust(wspace=0.3)
            axes = axes.flatten()
            for k, pk_key in enumerate(parameters):
                ax = axes[k]
                ax.axvline(x=1.0, color="grey", linestyle="--")

                # vertical lines
                if scan_key == "renal":
                    renal_colors = list(self.renal_colors.values())
                    for k, value in enumerate(self.renal_map.values()):
                        ax.axvline(x=value, color=renal_colors[k], linestyle="--")

                for renal_key in self.renal_map:
                    for cirrhosis_key in self.cirrhosis_map:
                        if (cirrhosis_key != list(self.cirrhosis_map.keys())[0]) and (renal_key != list(self.renal_map.keys())[0]):
                            # not interested in double changes
                            continue

                        sim_key = f"scan_{scan_key}_{cirrhosis_key}_{renal_key}"
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

                        if cirrhosis_key != list(self.cirrhosis_map.keys())[0]:
                            ax.plot(
                                x,
                                y,
                                marker="o",
                                linestyle="-",
                                color=self.cirrhosis_colors[cirrhosis_key],
                                markeredgecolor="black",
                                label=cirrhosis_key,
                            )
                        if renal_key != list(self.renal_map.keys())[0]:
                            if scan_key != "renal":
                                ax.plot(
                                    x,
                                    y,
                                    marker="o",
                                    linestyle="-",
                                    color=self.renal_colors[renal_key],
                                    markeredgecolor="black",
                                    label=renal_key,
                                )
                        if cirrhosis_key == list(self.cirrhosis_map.keys())[0] and renal_key == list(self.renal_map.keys())[0]:
                            ax.plot(
                                x,
                                y,
                                marker="o",
                                linestyle="-",
                                color="black",
                                markeredgecolor="black",
                                label="Healthy",
                            )

                ax.tick_params(axis="both", labelsize=SorafenibSimulationExperiment.tick_font_size)
                # ax.tick_params(axis="y", labelsize=SiSorafenibSimulationExperiment.tick_font_size)
                ax.set_xlabel(scan_data["label"], fontdict=SorafenibSimulationExperiment.font)
                ax.set_ylabel(
                    f"{self.pk_labels[pk_key]} [{self.pk_units[pk_key]}]",
                    fontdict=self.font,
                )
                ax.set_ylim(bottom=0.0)

                ax.set_xscale("log")
                ax.legend(fontsize=SorafenibSimulationExperiment.legend_font_size)
                ax.legend()

            figures[f"fig_{scan_key}"] = f
        return figures


if __name__ == "__main__":
    run_experiments(ParametersScan, output_dir="scan_parameters")
