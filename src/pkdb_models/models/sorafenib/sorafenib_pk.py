"""Sorafenib pharmacokinetics."""
import pandas as pd

from pkdb_analysis.pk.pharmacokinetics import TimecoursePK
from sbmlsim.result import XResult
from sbmlutils.console import console
from sbmlutils.log import get_logger

logger = get_logger(__name__)


def calculate_sorafenib_pk(
    experiment: "SorafenibSimulationExperiment",
    xres: XResult,
    tstart: float,
    tend: float,
) -> pd.DataFrame:
    """
    Calculate sorafenib parameters.
    """
    if (tstart is not None) and (tend is not None):
        # filter the time vector with tstart and tend
        da = xres.sel(_time=slice(tstart, tend))
        da = da.assign_coords(_time=da._time-tstart)
        da["time"] = da["time"] - tstart
        # console.print(da)

        xdata = XResult(
            xdataset=da,
            uinfo=xres.uinfo
        )
    else:
        xdata = xres

    # scanned dimension
    scandim = xdata._redop_dims()[0]
    Q_ = experiment.Q_

    # reads the initial dose from the results
    print(xdata)
    print(xdata.keys())
    dose_vec = Q_(xdata["PODOSE_sor"].values[0], xdata.uinfo["PODOSE_sor"])

    # calculate all pharmacokinetics
    pk_dicts = list()
    for k_dose, dose in enumerate(dose_vec):

        t_vec: Q_ = xdata.dim_mean("time")
        t_vec = Q_(t_vec.magnitude, xdata.uinfo["time"])
        c_vec = Q_(
            xdata["[Cve_sor]"].sel({scandim: k_dose}).values, xdata.uinfo["[Cve_sor]"]
        )
        dose = dose_vec[k_dose]
        dose_mmole = dose / experiment.Mr.sor
        tcpk = TimecoursePK(
            time=t_vec,
            concentration=c_vec,
            substance="sorafenib",
            dose=dose_mmole,
            ureg=experiment.ureg,
        )
        pk_dicts.append(tcpk.pk.to_dict())

    return pd.DataFrame(pk_dicts)
