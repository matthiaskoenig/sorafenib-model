"""Sorafenib intestine model."""
import numpy as np

from pkdb_models.models.sorafenib.models import annotations, templates
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *


class U(templates.U):
    """UnitDefinitions"""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")
    min_per_hr = UnitDefinition("min_per_hr", "min/hr")


_m = Model(
    "sorafenib_intestine",
    name="Model for sorafenib elimination in the small intestine",
    notes="""
    # Model for sorafenib absorption from the small intestine

    Model for enterohepatic uptake of sorafenib via OATP2B1 (SLCO2B1) 

    """
    + templates.terms_of_use,
    units=U,
    model_units=templates.model_units,
    creators=templates.creators,
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
    ),
    Compartment(
        "Vgu",
        1.2825,  # 0.0171 [l/kg] * 75 kg
        name="intestine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["gu"],
    ),
    Compartment(
        "Vlumen",
        1.2825 * 0.9,  # 0.0171 [l/kg] * 75 kg * 0.9,
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["stomach"],
    ),
    Compartment(
        "Ventero",
        np.nan,
        name="intestinal lining (enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        spatialDimensions=2,
        annotations=annotations.compartments["apical"],
        notes="""
        Membrane with microvilli towards the lumen.
        """
    ),
    Compartment(
        "Vbaso",
        np.nan,
        name="basolateral membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        spatialDimensions=2,
        annotations=annotations.compartments["basolateral"],
        notes="""
        Membrane towards the capillaries (blood).
        """
    ),
]

_m.species = [
    Species(
        "sor_stomach",
        metaId="meta_sor_stomach",
        initialConcentration=0.0,
        name="sorafenib (stomach)",
        compartment="Vstomach",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sor"],
        port=True,
    ),
    Species(
        "sor_lumen",
        initialConcentration=0.0,
        name="sorafenib (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sor"],
        port=True,
    ),
    Species(
        "sor_ext",
        initialConcentration=0.0,
        name="sorafenib (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sor"],
        port=True,
    ),
    Species(
        "sor_feces",
        initialConcentration=0.0,
        name="sorafenib (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sor"],
        port=True,
    ),
    Species(
        "sg_lumen",
        initialConcentration=0.0,
        name="sorafenib N-oxide glucuronide (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
        port=True,
    ),
    Species(
        "sg_ext",
        initialConcentration=0.0,
        name="sorafenib N-oxide glucuronide (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
        port=True,
    ),
    Species(
        "sg_feces",
        initialConcentration=0.0,
        name="sorafenib N-oxide glucuronide (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        f"F_sor_abs",
        0.5,
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"fraction absorbed sorafenib",
        notes="""
        Fraction absorbed, i.e., only a fraction of the sorafenib in the intestinal lumen
        is absorbed. This parameter determines how much of the sorafenib is excreted in
        every round of the enterohepatic circulation.

        Also, this parameter is dependent on OATP2B1 function, taking into consideration there may be
        SNPs hindering/enhancing the enzyme's function in the active process of the uptake of sorafenib .
        If OATP2B1 function is hindered, there would be more sorafenib  excreted into faeces than absorbed and vice
        versa. (added: 16.03.2022)

        `F_sor_abs` of dose is absorbed. `(1-F_sor_abs)` is excreted in feces.
        """,
    ),
    # Parameter(
    #     f"F_sg_abs",
    #     0.5,
    #     U.dimensionless,
    #     constant=True,
    #     sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    #     name=f"fraction absorbed SG",
    #     notes="""
    # Fraction absorbed, i.e., only a fraction of the SG in the intestinal lumen
    # is absorbed. This parameter determines how much of the SG is excreted in
    # every round of the enterohepatic circulation.
    #
    # `F_sg_abs` of dose is absgbed. `(1-F_sg_abs)` is excreted in feces.
    # """,
    # ),
    Parameter(
        "SORABS_Vmax",
        0.02,
        unit=U.per_min,
        name="Vmax for sorafenib absorption",
        sboTerm=SBO.MAXIMAL_VELOCITY,
    ),
    # Parameter(
    #     "SGABS_Vmax",
    #     0.02,
    #     unit=U.per_min,
    #     name="Vmax for SG absorption",
    #     sboTerm=SBO.MAXIMAL_VELOCITY,
    # ),
]

_m.rules.extend([
    AssignmentRule(
        "absorption_sor",
        value="SORABS_Vmax * Vgu * sor_lumen",
        unit=U.mmole_per_min,
        name="absorption SOR",
    ),
    AssignmentRule(
        "absorption_sg",
        value="SORABS_Vmax * Vgu * sg_lumen",
        unit=U.mmole_per_min,
        name="absorption SG",
        notes="""
        Assumption: identical Vmax for SG as for SOR.
        """
    )
])

_m.reactions = [
    Reaction(
        "SORABS",
        name="absorption sorafenib (OATP2B1)",
        equation="sor_lumen -> sor_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
            Parameter(
                sid="f_OATP2B1",
                name="parameter for OATP2B1 activity",
                value=1.0,
                unit=U.dimensionless,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
                parameter controls genetic variants (1.0 reference function,
                < 1.0 reduced function, > 1.0 increased function)
                """,
            ),
        ],
        formula=("F_sor_abs * absorption_sor * f_OATP2B1", U.mmole_per_min),
        notes="""
        Import via OATP2B1 is irreversible at apical membrane of enterocytes.
        """,
    ),
    # fraction excreted (not available for absorption)
    Reaction(
        sid=f"SOREXC",
        name=f"excretion sorafenib  (feces)",
        compartment="Vlumen",
        equation=f"sor_lumen -> sor_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            f"(1 dimensionless - F_sor_abs) * absorption_sor",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        "SGABS",
        name="absorption SG",
        equation="sg_lumen -> sg_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
        ],
        formula=("F_sor_abs * absorption_sg", U.mmole_per_min),
        notes="""
        Assumption: same fraction absorbed for SG as SOR to reduce parameters
        """
    ),
    # fraction excreted (not available for absorption)
    Reaction(
        sid=f"SGEXC",
        name=f"excretion SG (feces)",
        compartment="Vlumen",
        equation=f"sg_lumen -> sg_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            f"(1 dimensionless - F_sor_abs) * absorption_sg",
            U.mmole_per_min,
        ),
        notes="""
        Assumption: same fraction absorbed for SG as SOR to reduce parameters
        """
    ),
]


_m.parameters.extend([
    Parameter(
        f"PODOSE_sor",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose sorafenib [mg]",
        port=True,
    ),
    Parameter(
        f"POSTOMACH_sor",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose in stomach sor [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_application_sor",
        1000,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"Ka [1/hr] application sor",
        notes="""Fast application to shift applied dose in the stomach."""
    ),
    Parameter(
        f"Ka_dis_sor",
        2.0,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"Ka_dis [1/hr] dissolution sorafenib",
        port=True
    ),
    Parameter(
        f"Mr_sor",
        464.826,
        U.g_per_mole,
        constant=True,
        name=f"Molecular weight sorafenib [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
        port=True,
    ),
])

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"application_sor",
            name=f"application sorafenib",
            formula=(
                f"Ka_application_sor/60 min_per_hr * PODOSE_sor/Mr_sor",
                U.mmole_per_min,
            ),
            equation=f" -> sor_stomach",
            compartment="Vgu",
        ),

        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"dissolution_sor",
            name=f"dissolution sorafenib",
            formula=(
                f"Ka_dis_sor/60 min_per_hr * POSTOMACH_sor/Mr_sor",
                U.mmole_per_min,
            ),
            equation=f"sor_stomach -> sor_lumen",
            compartment="Vgu",
        ),
    ]
)
_m.rate_rules.extend(
    [
        RateRule(f"PODOSE_sor", f"-application_sor * Mr_sor", U.mg_per_min),
        RateRule(f"POSTOMACH_sor", f"(application_sor - dissolution_sor) * Mr_sor", U.mg_per_min),
    ]
)


model_intestine = _m

if __name__ == "__main__":
    from pkdb_models.models.sorafenib import MODEL_BASE_PATH

    result = create_model(
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        model=model_intestine, sbml_level=3, sbml_version=2
    )
    visualize_sbml(result.sbml_path, delete_session=False)
