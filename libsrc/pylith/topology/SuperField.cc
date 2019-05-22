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

#include <portinfo>

#include "SuperField.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*

#include <cassert> // USES assert()

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::topology::SuperField::SuperField(const pylith::topology::Field& field1,
                                         const pylith::topology::Field& field2) :
    _field1(field1),
    _field2(field2),
    _superIS(NULL),
    _superDM(NULL),
    _superVecLocal(NULL),
    _superVecGlobal(NULL),
    _field1VecGlobal(NULL),
    _field2VecGlobal(NULL) {}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::topology::SuperField::~SuperField(void) {
    deallocate();
} // destructor


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::SuperField::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;
    err = ISDestroy(&_superIS[0]);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&_superIS[1]);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&_superDM);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_superVecLocal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_superVecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_field1VecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecDestroy(&_field2VecGlobal);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ---------------------------------------------------------------------------------------------------------------------
// Get PETSc DM for super field.
PetscDM
pylith::topology::SuperField::getDM(void) const {
    return _superDM;
} // getDM


// ---------------------------------------------------------------------------------------------------------------------
// Get PETSc local vector for super field.
PetscVec
pylith::topology::SuperField::getLocalVector(void) const {
    return _superVecLocal;
} // getLocalVector


// ---------------------------------------------------------------------------------------------------------------------
// Initialize layout of super field.
void
pylith::topology::SuperField::initialize(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;

    PetscDM field1DM = _field1.dmMesh();
    PetscDM field2DM = _field2.dmMesh();

    // Create super field of [auxiliary_field, solution]
    PetscDM dms[2];
    dms[0] = field1DM;
    dms[1] = field2DM;
    err = DMCreateSuperDM(dms, 2, &_superIS, &_superDM);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(_superDM, &_superVecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMCreateLocalVector(_superDM, &_superVecLocal);PYLITH_CHECK_ERROR(err);

    // Create global vectors for fields.
    err = DMCreateGlobalVector(field1DM, &_field1VecGlobal);
    err = DMCreateGlobalVector(field2DM, &_field2VecGlobal);

    PYLITH_METHOD_END;
} // initialize


// ---------------------------------------------------------------------------------------------------------------------
// Copy values into super field.
void
pylith::topology::SuperField::copyTo(void) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = 0;

    PetscDM field1DM = _field1.dmMesh();
    PetscDM field2DM = _field2.dmMesh();

    // Copy field 1 to global vector
    err = DMLocalToGlobalBegin(field1DM, _field1.localVector(), INSERT_VALUES, _field1VecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(field1DM, _field1.localVector(), INSERT_VALUES, _field1VecGlobal);PYLITH_CHECK_ERROR(err);

    // Copy field 2 to global vector
    err = DMLocalToGlobalBegin(field2DM, _field2.localVector(), INSERT_VALUES, _field2VecGlobal);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(field2DM, _field2.localVector(), INSERT_VALUES, _field2VecGlobal);PYLITH_CHECK_ERROR(err);

    // Copy fields 1 and 2 to global super vector
    err = VecISCopy(_superVecGlobal, _superIS[0], SCATTER_FORWARD, _field1VecGlobal);PYLITH_CHECK_ERROR(err);
    err = VecISCopy(_superVecGlobal, _superIS[1], SCATTER_FORWARD, _field2VecGlobal);PYLITH_CHECK_ERROR(err);

    // Copy from global super vector to local super vector
    err = DMGlobalToLocalBegin(_superDM, _superVecGlobal, INSERT_VALUES, _superVecLocal);PYLITH_CHECK_ERROR(err);
    err = DMGlobalToLocalEnd(_superDM, _superVecGlobal, INSERT_VALUES, _superVecLocal);PYLITH_CHECK_ERROR(err);

    _field1.view("FIELD 1");
    _field2.view("FIELD 2");

    PYLITH_METHOD_END;
} // copyTo


// End of file
