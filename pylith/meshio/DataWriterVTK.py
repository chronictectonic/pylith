#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/DataWriterVTK.py
##
## @brief Python object for writing finite-element data to VTK file.
##
## Factory: output_data_writer

from DataWriter import DataWriter

# DataWriterVTK class
class DataWriterVTK(DataWriter):
  """
  Python object for writing finite-element data to VTK file.

  Factory: output_data_writer
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(DataWriter.Inventory):
    """
    Python object for managing DataWriterVTK facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing DataWriterVTK facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Name of VTK file.
    ## @li \b timeFormat C style format string for time stamp in filename.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="output.vtk")
    filename.meta['tip'] = "Name of VTK file."

    timeFormat = pyre.inventory.str("time_format", default="%f")
    timeFormat.meta['tip'] = "C style format string for time stamp in filename."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="solutioniovtk"):
    """
    Constructor.
    """
    DataWriter.__init__(self, name)
    return


  def initialize(self):
    """
    Initialize writer.
    """
    DataWriter.initialize(self)
    self.cppHandle.filename = self.filename
    self.cppHandle.timeFormat = self.timeFormat
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    DataWriter._configure(self)
    self.filename = self.inventory.filename
    self.timeFormat = self.inventory.timeFormat
    return


  def _createCppHandle(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.cppHandle:
      import pylith.meshio.meshio as bindings
      self.cppHandle = bindings.DataWriterVTK()
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def output_data_writer():
  """
  Factory associated with DataWriterVTK.
  """
  return DataWriterVTK()


# End of file 
