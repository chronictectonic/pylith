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

/** @file libsrc/bc/AbsorbingDampersNew.hh
 *
 * @brief C++ implementation of AbsorbingDampers for absorbing boundary conditions in dynamic problems.
 */

#if !defined(pylith_bc_absorbingDampersnew_hh)
#define pylith_bc_absorbingDampersnew_hh

// Include directives ---------------------------------------------------
#include "BoundaryConditionNew.hh" // ISA BoundaryCondition
#include "pylith/feassemble/IntegratorPointwise.hh" // ISA IntegratorPointwise

#include "pylith/topology/topologyfwd.hh" // USES Field

// AbsorbingDampersNew ----------------------------------------------------
/// @brief AbsorbingDampers (e.g., traction) boundary conditions.
class pylith::bc::AbsorbingDampersNew :
    public BoundaryConditionNew,
    public pylith::feassemble::IntegratorPointwise { // class AbsorbingDampersNew
    friend class TestAbsorbingDampersNew;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    AbsorbingDampersNew(void);

    /// Destructor.
    ~AbsorbingDampersNew(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void refDir1(const double vec[3]);

    /** Set second choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
     */
    void refDir2(const double vec[3]);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Initialize boundary condition.
     *
     * @param[in] solution Solution field.
     */
    void initialize(const pylith::topology::Field& solution);

    /** Compute RHS residual for G(t,s).
     *
     * @param[out] residual Field for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void computeRHSResidual(pylith::topology::Field* residual,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution);

    /** Compute RHS Jacobian and preconditioner for G(t,s).
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void computeRHSJacobian(PetscMat jacobianMat,
                            PetscMat preconMat,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution);

    /** Compute LHS residual for F(t,s,\dot{s}).
     *
     * @param[out] residual Field for residual.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    void computeLHSResidual(pylith::topology::Field* residual,
                            const PylithReal t,
                            const PylithReal dt,
                            const pylith::topology::Field& solution,
                            const pylith::topology::Field& solutionDot);

    /** Compute LHS Jacobian and preconditioner for F(t,s,\dot{s}) with implicit time-stepping.
     *
     * @param[out] jacobianMat PETSc Mat with Jacobian sparse matrix.
     * @param[out] precondMat PETSc Mat with Jacobian preconditioning sparse matrix.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] tshift Scale for time derivative.
     * @param[in] solution Field with current trial solution.
     * @param[in] solutionDot Field with time derivative of current trial solution.
     */
    void computeLHSJacobianImplicit(PetscMat jacobianMat,
                                    PetscMat precondMat,
                                    const PylithReal t,
                                    const PylithReal dt,
                                    const PylithReal tshift,
                                    const pylith::topology::Field& solution,
                                    const pylith::topology::Field& solutionDot);


    /** Compute inverse of lumped LHS Jacobian for F(t,s,\dot{s}) with explicit time-stepping.
     *
     * @param[out] jacobianInv Inverse of lumped Jacobian as a field.
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     * @param[in] solution Field with current trial solution.
     */
    void computeLHSJacobianLumpedInv(pylith::topology::Field* jacobianInv,
                                     const PylithReal t,
                                     const PylithReal dt,
                                     const pylith::topology::Field& solution);



    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void _auxFieldSetup(const pylith::topology::Field& solution);

    /** Set kernels for RHS residual G(t,s).
     *
     * Potentially, there are g0 and g1 kernels for each equation. If no
     * kernel is needed, then set the kernel function to NULL.
     *
     * @param solution Solution field.
     */
    virtual
    void _setFEKernelsRHSResidual(const pylith::topology::Field& solution) const;

    /** Set constants used in finite-element integrations.
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     */
    virtual
    void _setFEConstants(const pylith::topology::Field& solution,
                         const PylithReal dt) const;

    /** Get factory for setting up auxliary fields.
     *
     * @returns Factor for auxiliary fields.
     */
    pylith::feassemble::AuxiliaryFactory* _auxFactory(void);


    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::topology::Mesh* _boundaryMesh;   ///< Boundary mesh.
    pylith::topology::FieldBase::Description _description; ///< Description of field associated with BC.
    PylithReal _refDir1[3]; ///< First choice reference direction used to compute boundary tangential directions.
    PylithReal _refDir2[3]; ///< Second choice reference direction used to compute boundary tangential directions.

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    pylith::bc::AbsorbingDampersAuxiliaryFactory* _auxAbsorbingDampersFactory; ///< Factory for auxiliary subfields.
    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    AbsorbingDampersNew(const AbsorbingDampersNew&); ///< Not implemented.
    const AbsorbingDampersNew& operator=(const AbsorbingDampersNew&); ///< Not implemented.

}; // class AbsorbingDampersNew

#endif // pylith_bc_absorbingDampersnew_hh


// End of file