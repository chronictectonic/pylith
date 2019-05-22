// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/topology/SuperField.hh
 *
 * @brief C++ object for combining two fields into a single (super) field.
 */

#if !defined(pylith_topology_updatederivedfield_hh)
#define pylith_topology_updatederivedfield_hh

#include "pylith/topology/topologyfwd.hh"

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/petscfwd.h" // USES PetscIS, PetscDM, PetscVec

class pylith::topology::SuperField : public pylith::utils::GenericComponent {
    friend class TestSuperField; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor.
     *
     * @param[in] field1 Field 1.
     * @param[in] field2 Field 2.
     */
    SuperField(const pylith::topology::Field& field1,
               const pylith::topology::Field& field2);

    /// Destructor.
    virtual ~SuperField(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get PETSc DM for super field.
     *
     * @returns PETSc DM for super field.
     */
    PetscDM getDM(void) const;

    /** Get PETSc local vector for super field.
     *
     * @returns PETSc local vector for super field.
     */
    PetscVec getLocalVector(void) const;

    /// Initialize layout of super field.
    void initialize(void);

    /// Copy values into super field.
    void copyTo(void);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    const pylith::topology::Field& _field1; ///< First field in super field.
    const pylith::topology::Field& _field2; ///< Second field in super field.
    PetscIS* _superIS; ///< Petsc IS for mapping b/t auxiliary DM + field2 DM and super DM.
    PetscDM _superDM; ///< Petsc DM for super field.
    PetscVec _superVecLocal; ///< Petsc Vec with local vector for super field.
    PetscVec _superVecGlobal; ///< Petsc Vec with global vector for super field.
    PetscVec _field1VecGlobal; ///< Petsc Vec with global vector for field 1.
    PetscVec _field2VecGlobal; ///< Petsc Vec with global vector for field 2.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    SuperField(const SuperField &); ///< Not implemented
    const SuperField& operator=(const SuperField&); ///< Not implemented

}; // class SuperField

#endif // pylith_topology_updatederivedfield_hh

// End of file
