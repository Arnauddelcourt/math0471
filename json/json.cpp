// JSON simple example
// This example does not handle errors.

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
//#include <rapidjson/istreamwrapper.h>
#include "rapidjson/stringbuffer.h"
#include <iostream>
#include "rapidjson/filereadstream.h"
#include <cstdio>

using namespace rapidjson;

int main(int argc, char **argv) 
{
    if(argc!=2)
    {
        std::cout << "usage: " << argv[0] << " file.json" << std::endl;
        return 1;
    }
/*
    // not available in ubuntu 14.04LTS
    std::ifstream ifs("test1.json");
    IStreamWrapper isw(ifs);
    Document d;
    d.ParseStream(isw);
*/
    FILE* fp = fopen(argv[1], "rb"); // non-Windows use "r"
    if(!fp)
    {
        std::cout << "file not found!"<< std::endl;
        return 1;
    }
    char readBuffer[65536];
    FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    Document d;
    d.ParseStream(is);

    // stringify the DOM
    StringBuffer buffer;
    Writer<StringBuffer> writer(buffer);
    d.Accept(writer);
    // print DOM
    std::cout << buffer.GetString() << std::endl;
    
    return 0;
}