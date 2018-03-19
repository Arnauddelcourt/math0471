#include "json.h"

#include "rapidjson/document.h"
#include "rapidjson/error/en.h"
#include "rapidjson/writer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"
#include <iostream>
#include <cstdio>

JSON_API void read_json(std::string const &fname, rapidjson::Document &d)
{
    FILE* fp = fopen(fname.c_str(), "rb"); // non-Windows use "r"
    if(!fp)
    {
        std::cout << "file " << fname << " not found!"<< std::endl;
        abort();
    }

    fseek(fp, 0, SEEK_END);
    size_t length = static_cast<size_t>(ftell(fp));
    //std::cout << "file size = " << length << std::endl;
    fseek(fp, 0, SEEK_SET);
    char* readBuffer = static_cast<char*>(malloc(length + 1));
    size_t readLength = fread(readBuffer, 1, length, fp);
    readBuffer[readLength] = '\0';
    fclose(fp);

    d.Parse<rapidjson::kParseCommentsFlag>(readBuffer);

    if(d.HasParseError())
    {
        std::cout << "ERROR: json document cannot be parsed!" << std::endl;
        std::cout << "\t Parse Error " << d.GetParseError() << " (" << GetParseError_En(d.GetParseError()) << ")" << std::endl;
        std::cout << "\t Error offset = " << d.GetErrorOffset() << std::endl;
        abort();
    }    

    if(!d.IsObject())
    {
        std::cout << "ERROR: json document is not an object!" << std::endl;
        abort();
    }
}

JSON_API bool read_bool(rapidjson::Document const &d, char const *name, bool def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsBool())
        {
            std::cout << "ERROR: \"" << name << "\" should be a bool" << std::endl;
            abort(); 
        }
        return d[name].GetBool();
    }
    return def;
}

JSON_API int read_int(rapidjson::Document const &d, char const *name, int def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsInt())
        {
            std::cout << "ERROR: \"" << name << "\" should be an integer" << std::endl;
            abort(); 
        }
        return d[name].GetInt();
    }
    return def;
}

JSON_API double read_double(rapidjson::Document const &d, char const *name, double def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsNumber())
        {
            std::cout << "ERROR: \"" << name << "\" should be a number" << std::endl;
            abort(); 
        }
        return d[name].GetDouble();
    }
    return def;
}

JSON_API Vec3d read_Vec3d(rapidjson::Document const &d, char const *name, Vec3d const &def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsArray())
        {
            std::cout << "ERROR: \"" << name << "\" should be an array" << std::endl;
            abort(); 
        }
        if(d[name].Size()!=3)
        {
            std::cout << "ERROR: wrong array size for \"" << name <<"\"!" << std::endl;
            abort(); 
        }
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
        {
            if(!d[name][i].IsNumber())
            {
                std::cout << "ERROR: array \"" << name <<"\" does not contain numbers!" << std::endl;
                abort();
            }
        }
        return Vec3d(d[name][0].GetDouble(), d[name][1].GetDouble(), d[name][2].GetDouble());
    }
    return def;
}

JSON_API Vec3i read_Vec3i(rapidjson::Document const &d, char const *name, Vec3i const &def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsArray())
        {
            std::cout << "ERROR: \"" << name << "\" should be an array" << std::endl;
            abort(); 
        }
        if(d[name].Size()!=3)
        {
            std::cout << "ERROR: wrong array size for \"" << name <<"\"!" << std::endl;
            abort(); 
        }
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
        {
            if(!d[name][i].IsInt())
            {
                std::cout << "ERROR: array \"" << name <<"\" does not contain integers!" << std::endl;
                abort();
            }
        }
        return Vec3i(d[name][0].GetInt(), d[name][1].GetInt(), d[name][2].GetInt());
    }
    return def;
}
