"""FitParameters for parameter fitting."""

from sbmlsim.fit import FitParameter

parameters = [

    # dissolution in stomach sorafenib
    FitParameter(
        pid="Ka_dis_sor",
        start_value=2,
        lower_bound=0.1,
        upper_bound=10,
        unit="1/hr",
    ),
    # fraction absorbed in stomach SOR and SG
    FitParameter(
        pid="GU__F_sor_abs",
        start_value=0.5,
        lower_bound=0.4,
        upper_bound=0.6,
        unit="dimensionless",
    ),

    # absorption intestine (OATP2B1)
    FitParameter(
        pid="GU__SORABS_Vmax",
        start_value=0.2,
        lower_bound=1E-4,
        upper_bound=1E2,
        unit="1/min",
    ),

    # kidney excretion
    FitParameter(
        pid="KI__SGEX_k",
        start_value=0.01,
        lower_bound=1E-3,
        upper_bound=1,
        unit="1/min",
    ),


    # liver uptake
    FitParameter(
        pid="LI__SORIM_Vmax",
        start_value=1,
        lower_bound=1E-2,
        upper_bound=10,
        unit="mmole/min/l",
    ),
    FitParameter(
        pid="LI__M2EX_Vmax",
        start_value=1,
        lower_bound=1E-2,
        upper_bound=10,
        unit="mmole/min/l",
    ),
    FitParameter(
        pid="LI__SGEX_Vmax",
        start_value=1,
        lower_bound=1E-2,
        upper_bound=10,
        unit="mmole/min/l",
    ),
    FitParameter(
        pid="LI__SOR2M2_Vmax",
        start_value=0.1,
        lower_bound=1E-2,
        upper_bound=1,
        unit="mmole/min/l",
    ),
    FitParameter(
        pid="LI__M2GLU_Vmax",
        start_value=0.1,
        lower_bound=1E-2,
        upper_bound=1,
        unit="mmole/min/l",
    ),
]

parameters_all = parameters
