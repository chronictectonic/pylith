// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputSolnPoints.hh
 *
 * @brief C++ object for managing output of finite-element data over a
 * subdomain.
 */

#if !defined(pylith_meshio_outputsolnpoints_hh)
#define pylith_meshio_outputsolnpoints_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations
#include "pylith/utils/petscfwd.h"

#include "pylith/topology/Mesh.hh" // ISA OutputManager<Mesh>
#include "pylith/topology/Field.hh" // ISA OutputManager<Field<Mesh>>
#include "OutputManager.hh" // ISA OutputManager

// OutputSolnPoints -----------------------------------------------------
/** @brief C++ object for managing output of finite-element data over
 * a subdomain.
 */
class pylith::meshio::OutputSolnPoints : 
  public OutputManager<topology::Mesh, topology::Field<topology::Mesh> >
{ // OutputSolnPoints
  friend class TestOutputSolnPoints; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  OutputSolnPoints(void);

  /// Destructor
  ~OutputSolnPoints(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Get mesh associated with points.
   *
   * @returns Mesh associated with points.
   */
  const topology::Mesh& pointsMesh(void);

  /** Setup interpolator.
   *
   * @param mesh Domain mesh.
   * @param points Array of coordinates for points [numPoints*spaceDim].
   * @param numPoints Number of points.
   * @param spaceDim Spatial dimension for coordinates.
   */
  void setupInterpolator(topology::Mesh* mesh,
			 const PylithScalar* points,
			 const int numPoints,
			 const int spaceDim);
  
  /** Append finite-element vertex field to file.
   *
   * @param t Time associated with field.
   * @param field Vertex field.
   * @param mesh Mesh for output.
   */
  void appendVertexField(const PylithScalar t,
			 topology::Field<topology::Mesh>& field,
			 const topology::Mesh& mesh);

  /** Append finite-element cell field to file.
   *
   * @param t Time associated with field.
   * @param field Cell field.
   * @param label Name of label defining cells to include in output
   *   (=0 means use all cells in mesh).
   * @param labelId Value of label defining which cells to include.
   */
  void appendCellField(const PylithScalar t,
		       topology::Field<topology::Mesh>& field,
		       const char* label =0,
		       const int labelId =0);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  OutputSolnPoints(const OutputSolnPoints&); ///< Not implemented.
  const OutputSolnPoints& operator=(const OutputSolnPoints&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  topology::Mesh* _mesh; ///< Domain mesh.
  topology::Mesh* _pointsMesh; ///< Mesh for points (no cells).
  PetscDMMeshInterpolationInfo _interpolator; ///< Field interpolator.

}; // OutputSolnPoints

#endif // pylith_meshio_outputsolnpoints_hh

// End of file 