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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/** @file modulesrc/feassemble/GeometryQuad3D.i
 *
 * @brief Python interface to C++ GeometryQuad3D object.
 */

namespace pylith {
  namespace feassemble {

    class GeometryQuad3D : public CellGeometry
    { // GeometryQuad3D

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
  
      /// Default constructor.
      GeometryQuad3D(void);
      
      /// Default destructor.
      ~GeometryQuad3D(void);
      
      /** Create a copy of geometry.
       *
       * @returns Copy of geometry.
       */
      CellGeometry* clone(void) const;
      
      /** Get cell geometry for lower dimension cell.
       *
       * @returns Pointer to cell geometry object corresponding to next
       * lower dimension, NULL if there is no lower dimension object.
       */
      CellGeometry* geometryLowerDim(void) const;
      
      /** Transform coordinates in reference cell to global coordinates.
       *
       * @param ptsGlobal Array of points in global coordinate system.
       * @param ptsRef Array of points in reference cell.
       * @param vertices Array of cell vertices in global coordinates.
       * @param dim Dimension of global coordinate system.
       * @param npts Number of points to transform.
       */
      void ptsRefToGlobal(PylithScalar* ptsGlobal,
			  const PylithScalar* ptsRef,
			  const PylithScalar* vertices,
			  const int dim,
			  const int npts =1) const;
      
      /** Compute Jacobian at location in cell.
       *
       * @param jacobian Jacobian at location.
       * @param det Determinant of Jacobian at location.
       * @param vertices Coordinates of vertices of cell.
       * @param numVertices Number of vertices in cell.
       * @param spaceDim Spatial dimension of coordinates.
       * @param location Location in reference cell at which to compute Jacobian.
       * @param cellDim Dimension of reference cell.
       */
      void jacobian(pylith::scalar_array* jacobian,
		    PylithScalar* det,
		    const PylithScalar* vertices,
		    const int numVertices,
		    const int spaceDim,
		    const PylithScalar* location,
		    const int cellDim) const;

      /** Compute Jacobian at location in cell.
       *
       * @param jacobian Jacobian at location.
       * @param det Determinant of Jacobian at location.
       * @param vertices Coordinates of vertices of cell.
       * @param ptsRef Points in reference cell at which to compute Jacobian.
       * @param dim Dimension of coordinate system.
       * @param npts Number of points to transform.
       */
      void jacobian(PylithScalar* jacobian,
		    PylithScalar* det,
		    const PylithScalar* vertices,
		    const PylithScalar* ptsRef,
		    const int dim,
		    const int npts =1) const;
      
    }; // GeometryQuad3D

  } // feassemble
} // pylith

// End of file
