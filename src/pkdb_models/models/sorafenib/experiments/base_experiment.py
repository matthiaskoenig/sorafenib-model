"""Reusable functionality for multiple simulation experiments."""
from collections import namedtuple
from typing import Dict

import pandas as pd
import numpy as np
import matplotlib
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

from pkdb_models.models.sorafenib import MODEL_PATH

# Constants for conversion
from pkdb_models.models.sorafenib.sorafenib_pk import calculate_sorafenib_pk

MolecularWeights = namedtuple("MolecularWeights", "sor m2 sg")



class SorafenibSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 13
    suptitle_font_size = 40

    units: Dict[str, str] = {
        "time": "hr",

        "[Cve_sor]": "µM",
        "[Cve_m2]": "µM",
        "[Cve_sg]": "µM",

        "Aurine_sg": "µmole",
        "Afeces_sor": "µmole",
        "Afeces_sg": "µmole",

        "Cve_sg_sor": "dimensionless",
        "Cve_m2_sor": "dimensionless",
    }
    labels: Dict[str, str] = {
        "time": "time",

        "[Cve_sor]": "sorafenib (plasma)",
        "[Cve_m2]": "M2 (plasma)",
        "[Cve_sg]": "SG (plasma)",

        "Aurine_sg": "SG urine",
        "Afeces_sor": "sorafenib feces",
        "Afeces_sg": "SG feces",

        "Cve_sg_sor": "SG / sorafenib",
        "Cve_m2_sor": "M2 / sorafenib",
    }

    label_time = "time"
    label_sor = "sorafenib (plasma)"
    label_m2 = "M2 (plasma)"
    label_sg = "SG (plasma)"
    label_sg_urine = "SG urine"
    label_sor_feces = "sorafenib feces"
    label_sg_feces = "SG feces"

    pk_labels = {
        "auc": "AUCend",
        "aucinf": "AUC",
        "cl": "Total clearance",
        # "cl_renal": "Renal clearance",
        # "cl_hepatic": "Hepatic clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "tmax": "tmax",
        "kel": "kel",
    }

    pk_units = {
        "auc": "mmole/l*hr", #mmole/l⋅hr
        "aucinf": "mmole/l*hr",
        "cl": "l/min",
        # "cl_renal": "l/min",
        # "cl_hepatic": "l/min",
        "cmax": "µmole/l",
        "thalf": "hr",
        "tmax": "hr",
        "kel": "1/hr",
    }

    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,  # CPT-A
        "Moderate cirrhosis": 0.6979591836734694,  # CPT-B
        "Severe cirrhosis": 0.8127551020408164,  # CPT-C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#74a9cf",
        "Moderate cirrhosis": "#2b8cbe",
        "Severe cirrhosis": "#045a8d",
    }

    renal_map = {
        "Normal renal function": 101.0 / 101.0,  # 1.0,
        "Mild renal impairment": 69.5 / 101.0,  # 0.69
        "Moderate renal impairment": 32.5 / 101.0,  # 0.32
        "Severe renal impairment": 19.5 / 101.0,  # 0.19
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {
            # first optimization with Huang2017
            'Ka_dis_sor': Q_(0.6380404886582229, '1/hr'),  # [0.1 - 10]
            'GU__F_sor_abs': Q_(0.5999999998431342, 'dimensionless'),  # [0.4 - 0.6]
            'GU__SORABS_Vmax': Q_(0.0004564619167017758, '1/min'),  # [0.0001 - 100.0]
            'KI__SGEX_k': Q_(0.00100102340272109, '1/min'),  # [0.001 - 1]
            'LI__SORIM_Vmax': Q_(0.01000000148281594, 'mmole/min/l'),  # [0.01 - 10]
            'LI__M2EX_Vmax': Q_(0.043035916104154845, 'mmole/min/l'),  # [0.01 - 10]
            'LI__SGEX_Vmax': Q_(4.852607840633812, 'mmole/min/l'),  # [0.01 - 10]
            'LI__SOR2M2_Vmax': Q_(0.010000000008669747, 'mmole/min/l'),  # [0.01 - 1]
            'LI__M2GLU_Vmax': Q_(0.021752867459769464, 'mmole/min/l'),  # [0.01 - 1]

        }
        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return SorafenibSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",
                "[Cve_sor]",
                "[Cve_m2]",
                "[Cve_sg]",
                "Afeces_sor",
                "Aurine_sg",
                "Afeces_sg",
                "Cve_sg_sor",
                "Cve_m2_sor",
                "PODOSE_sor",
                "f_cirrhosis",
                "KI__f_renal_function",
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            sor=self.Q_(464.826, "g/mole"),  # 50 mg/464 g/mole
            m2=self.Q_(421.6, "g/mole"),
            sg=self.Q_(472.6, "g/mole"),
        )

    def calculate_sorafenib_pk(self, scans: list = [], tstart: float = None, tend: float = None) -> Dict[str, pd.DataFrame]:
        """Calculate sorafenib parameters for simulations (scans)"""
        pk_dfs = {}
        if scans:
            for sim_key in scans:
                xres = self.results[f"task_{sim_key}"]
                df = calculate_sorafenib_pk(experiment=self, xres=xres, tstart=tstart, tend=tend)
                pk_dfs[sim_key] = df
        else:
            for sim_key in self._simulations.keys():
                xres = self.results[f"task_{sim_key}"]
                df = calculate_sorafenib_pk(experiment=self, xres=xres, tstart=tstart, tend=tend)
                pk_dfs[sim_key] = df
        return pk_dfs

    @staticmethod
    def color_for_dose(dose: float):
        """Colors for given sorafenib dose."""

        dose_max: float = 2800.0
        # dose_relative = np.log(dose) / np.log(dose_max)
        dose_relative: float = dose / dose_max

        cmap = matplotlib.cm.get_cmap('plasma_r')
        # cmap = matplotlib.cm.get_cmap('inferno_r')
        if np.isclose(dose, 0.0):
            # use black for placebo
            color_rgba = "#000000"
        else:
            color_rgba = cmap(dose_relative)
        return color_rgba
