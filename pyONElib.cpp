#include <cstring>
#include "ONElib.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

// Python passes line types as one-character strings
static char lineTypeChar (const std::string &lineType)
{
    if (lineType.length() != 1) throw std::invalid_argument("lineType must be a single character");
    return lineType[0];
}

PYBIND11_MODULE(ONEcode, m) {
    py::class_<ONEschema>(m, "ONEschema")
        .def(py::init<const std::string &>());

    py::class_<ONEfile>(m, "ONEfile")
        .def(py::init<const std::string &>())
        .def(py::init<const std::string &, int>())
        .def(py::init<const std::string &, const std::string &, const ONEschema &, const std::string &, int>())
        .def(py::init<const std::string &, const std::string &, ONEfile &, int>())
        .def("checkSchemaText", &ONEfile::checkSchemaText)
        // the type of the line read, or "" at the end of the file, so that "while f.readLine():" terminates
        .def("readLine", [](ONEfile &self) {
            char t = self.readLine();
            return t ? std::string(1, t) : std::string();
        })
        .def("length", &ONEfile::listLength)
        .def("getInt", &ONEfile::getInt)
        .def("setInt", &ONEfile::setInt)
        .def("getReal", &ONEfile::getReal)
        .def("setReal", &ONEfile::setReal)
        .def("getChar", &ONEfile::getChar)
        .def("setChar", &ONEfile::setChar)
        .def("getIntList", &ONEfile::getIntList)
        .def("getRealList", &ONEfile::getRealList)
        .def("getDNAchar", &ONEfile::getDNAchar)
        .def("getDNA2bit", &ONEfile::getDNA2bit)
        .def("getString", &ONEfile::getString)
        .def("getStringList", &ONEfile::getStringList)
        .def("getComment", &ONEfile::getComment)
        .def("gotoObject", &ONEfile::gotoObject)
        .def("lineType", &ONEfile::lineType)
        .def("lineNumber", &ONEfile::lineNumber)
        // number, maximum length and total length of lines of a type: as given in the header when
        // reading, and as accumulated so far when writing
        .def("count", [](ONEfile &self, const std::string &t) { return self.count(lineTypeChar(t)); })
        .def("max",   [](ONEfile &self, const std::string &t) { return self.max(lineTypeChar(t)); })
        .def("total", [](ONEfile &self, const std::string &t) { return self.total(lineTypeChar(t)); });
}
