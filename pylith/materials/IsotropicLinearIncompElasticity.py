# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/materials/IsotropicLinearIncompElasticity.py
#
# @brief Python material for isotropic, linearly elastic, incompressible
# material.
#
# Factory: incompressible_elasticity_rheology

from .RheologyIncompressibleElasticity import RheologyIncompressibleElasticity
from .materials import IsotropicLinearIncompElasticity as ModuleLinearElasticity


class IsotropicLinearIncompElasticity(RheologyIncompressibleElasticity, ModuleLinearElasticity):
    """
    Python material for isotropic, linearly elastic incompressible.

    INVENTORY

    Properties
      - *use_reference_state* Use reference stress/strain state.

    Facilities
      - *auxiliary_subfields* Discretization of physical properties and state variables.

    FACTORY: material
    """

    import pyre.inventory

    useReferenceState = pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    from .AuxSubfieldsIsotropicLinearElasticity import AuxSubfieldsIsotropicLinearElasticity
    from pylith.topology.Subfield import subfieldFactory
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxSubfieldsIsotropicLinearElasticity)
    auxiliarySubfields.meta['tip'] = "Discretization of physical properties and state variables."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearincompelasticity"):
        """
        Constructor.
        """
        RheologyIncompressibleElasticity.__init__(self, name)
        return

    def preinitialize(self, mesh):
        RheologyIncompressibleElasticity.preinitialize(self, mesh)

        print(self)
        ModuleLinearElasticity.useReferenceState(self, self.useReferenceState)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        ModuleLinearElasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def incompressible_elasticity_rheology():
    """
    Factory associated with IsotropicLinearIncompElasticity.
    """
    return IsotropicLinearIncompElasticity()


# End of file
