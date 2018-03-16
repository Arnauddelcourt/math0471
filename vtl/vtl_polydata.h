#ifndef VTL_POLYDATA_H
#define VTL_POLYDATA_H

#include "vtl.h"

void export_polydata_LEGACY(std::string const &filename,
                            int step,
                            std::vector<double> const &pos,
                            std::map<std::string, std::vector<double> *> const &scalars,
                            std::map<std::string, std::vector<double> *> const &vectors,
                            bool binary);

void export_polydata_XML(std::string const &filename,
                         int step,
                         std::vector<double> const &pos,
                         std::map<std::string, std::vector<double> *> const &scalars,
                         std::map<std::string, std::vector<double> *> const &vectors,
                         bool binary,
                         bool usez);

VTL_API void export_polydata(std::string const &filename,
                                  int step,
                                  std::vector<double> const &pos,
                                  std::map<std::string, std::vector<double> *> const &scalars,
                                  std::map<std::string, std::vector<double> *> const &vectors,
                                  PFormat format);

#endif // VTL_POLYDATA_H