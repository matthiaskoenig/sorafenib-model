from dataclasses import dataclass
from enum import Enum

from sbmlsim.fit.objects import MappingMetaData


class Tissue(str, Enum):
    PLASMA = "plasma"
    SERUM = "serum"
    URINE = "urine"
    FECES = "feces"
    BILE = "bile"


class Route(str, Enum):
    PO = "po"
    IV = "iv"
    SL = "sl"


class Dosing(str, Enum):
    SINGLE = "single"
    MULTIPLE = "multiple"


class Health(str, Enum):
    NR = "not reported"
    HEALTHY = "healthy"
    RENAL_IMPAIRMENT = "renal impairment"
    CARDIAC_IMPAIRMENT = "cardiac impairment"
    LIVER_IMPAIRMENT = "liver impairment"


class Fasting(str, Enum):
    NR = "not reported"
    FASTING = "fasting"
    NONFASTING = "nonfasting"
    LIGHT_BR = "light breakfast"

class PKPDData(str, Enum):
    PK = "pk"
    PD = "pd"


@dataclass
class SorafenibMappingMetaData(MappingMetaData):
    """Metadata for fitting experiment."""

    tissue: Tissue
    route: Route
    dosing: Dosing
    health: Health
    fasting: Fasting
    data: PKPDData
    outlier: bool = False

    def to_dict(self):
        """Convert to dictionary."""
        return {
            "tissue": self.tissue.name,
            "route": self.route.name,
            "dosing": self.dosing.name,
            "health": self.health.name,
            "fasting": self.fasting.name,
            "data": self.data.name,
            "outlier": self.outlier,
        }