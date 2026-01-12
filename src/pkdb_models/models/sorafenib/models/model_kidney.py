"""Sorafenib kidney model."""
from pkdb_models.models.sorafenib.models import annotations, templates
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *


class U(templates.U):
    """UnitDefinitions"""
    ml_per_min_bsa = UnitDefinition("ml_per_min_bsa", "ml/min/(1.73 * m^2)")

version = 1
_m = Model(
    sid="sorafenib_kidney",
    name="Model for sorafenib elimination in the kidneys",
    notes=f"""
    # Model for sorafenib elimination in the kidneys

    Model for renal excretion of sorafenib metabolites. 
    encoded in <a href="http://sbml.org">SBML</a> format.
    
    **version** {version}
    """
    + templates.terms_of_use,
    units=U,
    model_units=templates.model_units,
    creators=templates.creators,
)

_m.compartments = [
    Compartment(
        sid="Vext",
        value=1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
    ),
    Compartment(
        "Vki",
        0.3,
        name="kidney",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["ki"],
    ),
    Compartment(
        "Vurine",
        1.0,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["urine"],
    ),
]

_m.species = [
    Species(
        "sg_ext",
        initialConcentration=0.0,
        name="sorafenib glucuronide (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
        port=True,
    ),
    Species(
        "sg_urine",
        initialConcentration=0.0,
        name="sorafenib glucuronide (urine)",
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["sg"],
        port=True,
        notes="""
        Urinary species are in amounts.
        """,
    ),
]

_m.parameters = [
    Parameter(
        sid="f_renal_function",
        name="parameter for kidney function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Parameter controls the health state/renal function
        (1.0 reference function, < 1.0 reduced function, > 1.0 increased function)
        """,
    )
]

_m.reactions = [
    Reaction(
        "SGEX",
        name="excretion sorafenib glucuronide",
        equation="sg_ext -> sg_urine",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vki",
        pars=[
            Parameter(
                "SGEX_k",
                0.01,
                unit=U.per_min,
                name="rate urinary excretion sorafenib",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("f_renal_function * SGEX_k * Vki * sg_ext", U.mmole_per_min),
        annotations=[(BQB.IS, "ncit/C75913")],  # renal clearance
    ),
]

model_kidney = _m


if __name__ == "__main__":
    from pkdb_models.models.sorafenib import MODEL_BASE_PATH

    result = create_model(
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        model=model_kidney, sbml_level=3, sbml_version=2
    )
    visualize_sbml(result.sbml_path, delete_session=True)
