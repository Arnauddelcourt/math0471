#ifndef READPARS_H
#define READPARS_H

#include "vtl.h"
#include "vtlVec3.h"
#include "rapidjson/document.h"
#include <string>

void read_json(std::string const &fname, rapidjson::Document &d);
bool read_bool(rapidjson::Document const &d, char const *name, bool def);
Vec3d read_Vec3d(rapidjson::Document const &d, char const *name, Vec3d const &def);
Vec3i read_Vec3i(rapidjson::Document const &d, char const *name, Vec3i const &def);

#endif //READPARS_H