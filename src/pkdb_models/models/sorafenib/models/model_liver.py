"""Liver model of sorafenib.

Hepatic metabolism:
- oxidation to sorafenib N-oxide (M2) via CYP3A4
- glucuronidation to sorafenib N-oxide glucuronide (M2 glucuronide) via UGT1A9
- possible other minor pathways
- enterohepatic circulation of M2 glucuronide


M2 metabolite/Sorafenib N-oxide: Huang2017, Ferrario2016
Sorafenib-glucuronide (SG): Hussaarts2020, Inaba2019

TODO: add glucuronidation reation: M2 -> SG & transport of SG

"""

import numpy as np

from pkdb_models.models.sorafenib.models import annotations, templates
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *


class U(templates.U):
    """UnitDefinitions"""

    pass


version = 1
_m = Model(
    sid="sorafenib_liver",
    name="Model for hepatic sorafenib metabolism.",
    notes=f"""
    Model for sorafenib metabolism.
    **version** {version}

    ## Changelog

    **version 1**

    - initial model for demonstration

    **version 2**

    - added conversion and export of sorafenib metabolites

    """
    + templates.terms_of_use,
    units=U,
    model_units=templates.model_units,
    creators=templates.creators,
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
    Compartment(
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True
    ),
    Compartment(
        "Vbaso",
        np.nan,
        name="basolateral membrane",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["basolateral"],
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["apical"],
    ),
    Compartment(
        "Vbi",
        1.0,
        name="bile",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["bi"],
        port=True,
    ),
    Compartment(
        "Vlumen",
        0.1,
        name="intestinal lumen",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),

]

_m.species = [
    Species(
        "sor_ext",
        name="sorafenib(plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sor"],
        port=True
    ),
      Species(
        "m2_ext",
        name="sorafenib N-oxdide M2 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
        port=True
    ),
    # Species(
    #     "sg_ext",
    #     name="sorafenib N-oxdide glucuronide (plasma)",
    #     initialConcentration=0.0,
    #     compartment="Vext",
    #     substanceUnit=U.mmole,
    #     hasOnlySubstanceUnits=False,  # this is a concentration
    #     sboTerm=SBO.SIMPLE_CHEMICAL,
    #     annotations=annotations.species["sg"],
    #     port=True
    # ),
    Species(
        "sor",
        name="sorafenib (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sor"],
    ),
   Species(
        "m2",
        name="sorafenib N-oxdide M2 (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m2"],
    ),
    Species(
        "sg",
        name="sorafenib N-oxdide glucuronide (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
    ),
    Species(
        "sg_bi",
        initialConcentration=0.0,
        name="sorafenib N-oxdide glucuronide (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
        notes="""
        Bile SG in amount.
        """,
    ),
    Species(
        "sg_lumen",
        initialConcentration=0.0,
        name="sorafenib N-oxdide glucuronide (lumen)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
        port=True,
    )
]

_m.reactions = [
    Reaction(
        sid="SORIM",
        name="sorafenib import OATP1B/OCT1 (SORIM)",
        equation="sor_ext <-> sor",
        compartment="Vbaso",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "SORIM_Vmax",
                1.0,
                U.mmole_per_min_l,
                name="Vmax sorafenib import",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "SORIM_Km_sor",
                0.1,
                U.mM,
                name="Km for sorafenib import",
                sboTerm=SBO.MICHAELIS_CONSTANT,
            )
        ],
        formula=(
            "SORIM_Vmax/SORIM_Km_sor * Vli * (sor_ext - sor)/(1 dimensionless + sor_ext/SORIM_Km_sor + sor/SORIM_Km_sor)"
        )
    ),

    Reaction(
        sid="SOR2M2",
        name="sorafenib oxdidation CYP3A4 (SOR2M2)",
        equation="sor -> m2",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "SOR2M2_Vmax",
                0.1,
                U.mmole_per_min_l,
                name="Vmax sorafenib model conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "SOR2M2_Km_sor",
                0.1,
                U.mM,
                name="Km sorafenib conversion",
                sboTerm=SBO.MICHAELIS_CONSTANT,
            )
        ],
        formula=(
            "SOR2M2_Vmax/SOR2M2_Km_sor * Vli * sor/(1 dimensionless + sor/SOR2M2_Km_sor)"
        ),
        notes="""phase I oxidation mediated by cytochrome P450 3A4 (CYP3A4)
        Catalyzed partially via oxidases.
        """
    ),

    Reaction(
        sid="M2EX",
        name="M2 export (M2EX)",
        equation="m2 <-> m2_ext",
        compartment="Vapical",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M2EX_Vmax",
                1.0,
                U.mmole_per_min_l,
                name="Vmax sorafenib export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "M2EX_Km_m2",
                0.1,
                U.mM,
                name="Km sorafenib export",
                sboTerm=SBO.MICHAELIS_CONSTANT,
            )
        ],
        formula=(
            "M2EX_Vmax/M2EX_Km_m2 * Vli * (m2 - m2_ext)/(1 dimensionless + m2_ext/M2EX_Km_m2 + m2/M2EX_Km_m2)"
        )
    ),
    Reaction(
        sid="M2GLU",
        name="M2 glucuronidation UGT1A9 (M2GLU)",
        equation="m2 -> sg",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "M2GLU_Vmax",
                0.1,
                U.mmole_per_min_l,
                name="Vmax sorafenib model conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "M2GLU_Km_m2",
                0.1,
                U.mM,
                name="Km sorafenib conversion",
                sboTerm=SBO.MICHAELIS_CONSTANT,
            )
        ],
        formula=(
            "M2GLU_Vmax * Vli * m2/(m2 + M2GLU_Km_m2)"
        ),
        notes="""
        phase II conjugation mediated by UDP glucuronosyltransferase 1A9 (UGT1A9)
        """
    ),
    Reaction(
        sid="SGEX",
        name="sg export ABCC2/ABCG2 (SGEX)",
        equation="sg -> sg_bi",
        compartment="Vapical",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "SGEX_Vmax",
                1.0,
                U.mmole_per_min_l,
                name="Vmax sg export to bile",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "SGEX_Km_sg",
                0.1,
                U.mM,
                name="Km sg export to bile",
                sboTerm=SBO.MICHAELIS_CONSTANT,
            )
        ],
        formula=(
            "SGEX_Vmax * Vli * sg/(sg + SGEX_Km_sg)"
        )
    ),

    Reaction(
        "SGEHC",
        name="sg enterohepatic circulation",
        equation="sg_bi -> sg_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vlumen",
        notes="""
        Transport of sg_bi via bile to the intestinal lumen. This is the
        enterohepatic circulation of sg-bi 

        Irreversible transport of sg_bi via Mass-Action-Kinetics,
        """,
        pars=[],
        formula=("SGEX", U.mmole_per_min),
    ),
]

model_liver = _m


if __name__ == "__main__":
    from pkdb_models.models.gdeobdtpa import MODEL_BASE_PATH

    result = create_model(
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        model=model_liver, sbml_level=3, sbml_version=2
    )
    visualize_sbml(result.sbml_path, delete_session=True)

