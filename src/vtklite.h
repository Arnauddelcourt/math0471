#ifndef VTKLITE_H
#define VTKLITE_H

#include <string>
#include <vector>
#include <map>

enum PFormat
{
    LEGACY_TXT = 0,
    LEGACY_BIN = 1,
    XML_BIN = 2,
    XML_BINZ = 3
};

void export_polydata(std::string const &filename,
                     int step,
                     std::vector<double> const &pos,
                     std::map<std::string, std::vector<double> *> const &scalars,
                     std::map<std::string, std::vector<double> *> const &vectors,
                     PFormat format);

void export_spoints(std::string const &filename,
                    int step,
                    double o[3],
                    double dx[3],
                    int np[3],
                    std::map<std::string, std::vector<double> *> const &scalars,
                    std::map<std::string, std::vector<double> *> const &vectors,
                    PFormat format);

#endif // VTKLITE_H
