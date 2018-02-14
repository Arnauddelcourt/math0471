#ifndef VTKLITE_H
#define VTKLITE_H

#include <string>
#include <vector>
#include <map>

#if defined(WIN32)
#ifdef vtl_EXPORTS
#define VTL_API __declspec(dllexport)
#else
#define VTL_API __declspec(dllimport)
#endif
#else
#define VTL_API
#endif

namespace vtl 
{
    class SPoints;
    class UPoints;
}

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
                    PFormat format,
                    int myid = -1,
                    int *myn1 = NULL,
                    int *myn2 = NULL);

void export_spoints_XMLP(std::string const &filename,
                        int step,
                        double o[3],
                        double dx[3],
                        int np[3],
                        std::map<std::string, std::vector<double> *> const &scalars,
                        std::map<std::string, std::vector<double> *> const &vectors,
                        bool binary,
                        bool usez,
                        std::vector< std::vector<int> > const &extents
                        );

#endif // VTKLITE_H