/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FILEIO_VTK_PVDFILE_H
#define FILEIO_VTK_PVDFILE_H

#include <string>
#include <vector>

namespace FileIO
{

/*! Writes a basic PVD file for use with Paraview.
 *
 */
class PVDFile
{
public:
    //! Set a PVD file path
    explicit PVDFile(std::string const& pvd_fname) : _pvd_filename(pvd_fname) {}

    //! Add a VTU file to this PVD file.
    void addVTUFile(std::string const& vtu_fname, double timestep);

private:
    std::string const _pvd_filename;
    std::vector<std::pair<double, std::string>> _datasets; // a vector of (time, VTU file name)
};

} // namespace FileIO

#endif // FILEIO_VTK_PVDFILE_H
