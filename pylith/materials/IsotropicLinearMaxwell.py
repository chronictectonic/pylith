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
# @file pylith/materials/IsotropicLinearMaxwell.py
#
# @brief Python material for isotropic, linearly Maxwell viscoelastic material.
#
# Factory: elasticity_rheology

from .RheologyElasticity import RheologyElasticity
from .materials import IsotropicLinearMaxwell as ModuleLinearElasticity


class IsotropicLinearMaxwell(RheologyElasticity, ModuleLinearElasticity):
    """
    Python material for isotropic, linear Maxwell viscoelastic.

    INVENTORY

    Properties
      - *use_reference_state* Use reference stress/strain state.

    Facilities
      - *auxiliary_subfields* Discretization of physical properties and state variables.

    FACTORY: elasticity_rheology
    """

    import pyre.inventory

    useReferenceState = pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    from .AuxSubfieldsIsotropicLinearMaxwell import AuxSubfieldsIsotropicLinearMaxwell
    from pylith.topology.Subfield import subfieldFactory
    auxiliarySubfields = pyre.inventory.facilityArray(
        "auxiliary_subfields", itemFactory=subfieldFactory, factory=AuxSubfieldsIsotropicLinearMaxwell)
    auxiliarySubfields.meta['tip'] = "Discretization of physical properties and state variables."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearmaxwell"):
        """
        Constructor.
        """
        RheologyElasticity.__init__(self, name)
        return

    def preinitialize(self, mesh):
        RheologyElasticity.preinitialize(self, mesh)

        ModuleLinearElasticity.useReferenceState(self, self.useReferenceState)

        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """
        Call constructor for module object for access to C++ object.
        """
        ModuleLinearElasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def elasticity_rheology():
    """
    Factory associated with IsotropicLinearMaxwell.
    """
    return IsotropicLinearMaxwell()


# End of file
