// This line compiled successfully:
// g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) pyONElib.cpp -o ONEcode$(python3-config --extension-suffix) ONElib.o

#include <cstring>
#include "ONElib.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(ONEcode, m) {
    py::class_<ONEschema>(m, "ONEschema")
        .def(py::init<const std::string &>());

    py::class_<ONEfile>(m, "ONEfile")
        .def(py::init<const std::string &, const std::string &, const ONEschema &, const std::string &, int>())
        .def(py::init<const std::string &, const std::string &, ONEfile &, int>())
        .def("checkSchemaText", &ONEfile::checkSchemaText)
        .def("readLine", &ONEfile::readLine)
        .def("length", &ONEfile::length)
        .def("getInt", &ONEfile::getInt)
        .def("setInt", &ONEfile::setInt)
        .def("getReal", &ONEfile::getReal)
        .def("setReal", &ONEfile::setReal)
        .def("getChar", &ONEfile::getChar)
        .def("setChar", &ONEfile::setChar)
        .def("setDNAchar", &ONEfile::setDNAchar)
        .def("getIntList", &ONEfile::getIntList)
        .def("getRealList", &ONEfile::getRealList)
        .def("getDNA2bit", &ONEfile::getDNA2bit)
        .def("getString", &ONEfile::getString)
        .def("nextString", &ONEfile::nextString)
        .def("getComment", &ONEfile::getComment)
        .def("gotoObject", &ONEfile::gotoObject)
        .def("gotoGroup", &ONEfile::gotoGroup)
        .def("lineType", &ONEfile::lineType)
        .def("lineNumber", &ONEfile::lineNumber)
        .def("object", &ONEfile::object)
        .def("group", &ONEfile::group)
        .def("count", &ONEfile::count)
        .def("max", &ONEfile::max)
        .def("total", &ONEfile::total)
        .def("groupCount", &ONEfile::groupCount)
        .def("groupTotal", &ONEfile::groupTotal);
}

