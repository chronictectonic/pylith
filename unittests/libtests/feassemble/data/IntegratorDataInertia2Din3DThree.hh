// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application integratorinertia2din3dthree.

#if !defined(pylith_feassemble_integratordatainertia2din3dthree_hh)
#define pylith_feassemble_integratordatainertia2din3dthree_hh

#include "IntegratorData.hh"

namespace pylith {
  namespace feassemble {
     class IntegratorDataInertia2Din3DThree;
  } // pylith
} // feassemble

class pylith::feassemble::IntegratorDataInertia2Din3DThree : public IntegratorData
{

public: 

  /// Constructor
  IntegratorDataInertia2Din3DThree(void);

  /// Destructor
  ~IntegratorDataInertia2Din3DThree(void);

private:

  static const int _numVertices;

  static const int _spaceDim;

  static const int _numCells;

  static const int _cellDim;

  static const int _numBasis;

  static const int _numQuadPts;

  static const int _fiberDim;

  static const double _vertices[];

  static const int _cells[];

  static const double _quadPts[];

  static const double _quadWts[];

  static const double _basis[];

  static const double _basisDeriv[];

  static const double _fieldIn[];

  static const double _valsAction[];

  static const double _valsMatrix[];

  static const double _valsLumped[];

};

#endif // pylith_feassemble_integratordatainertia2din3dthree_hh

// End of file
