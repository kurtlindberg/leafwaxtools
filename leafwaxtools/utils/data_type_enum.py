"""
This module enumarates each plant wax data type (e.g., n-alkanoic acids,
n-alkanes) to be called in other leafwaxtools apis.
"""

from enum import Enum

class DataType(Enum):

    NALKANOIC_ACID = 1
    NALKANE = 2
    # NALCOHOL = 3
