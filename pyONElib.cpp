#include <cstring>
#include "ONElib.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(ONEcode, m) {
    py::class_<ONEschema>(m, "ONEschema")
        .def(py::init<const std::string &>());

    py::class_<ONEfile>(m, "ONEfile")
        .def(py::init<const std::string &>())
        .def(py::init<const std::string &, int>())
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
        .def("getIntList", &ONEfile::getIntList)
        .def("getRealList", &ONEfile::getRealList)
        .def("getDNAchar", &ONEfile::getDNAchar)
        .def("getDNA2bit", &ONEfile::getDNA2bit)
        .def("getString", &ONEfile::getString)
        .def("nextString", &ONEfile::nextString)
        .def("getComment", &ONEfile::getComment)
        .def("gotoObject", &ONEfile::gotoObject)
        .def("lineType", &ONEfile::lineType)
        .def("lineNumber", &ONEfile::lineNumber)
        
        // Updated function bindings to handle string-to-char conversion
        .def("givenCount", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.givenCount(static_cast<unsigned char>(lineType[0]));
        })
        
        .def("givenMax", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.givenMax(static_cast<unsigned char>(lineType[0]));
        })

        .def("givenTotal", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.givenTotal(static_cast<unsigned char>(lineType[0]));
        })

        .def("currentCount", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.currentCount(static_cast<unsigned char>(lineType[0]));
        })

        .def("currentMax", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.currentMax(static_cast<unsigned char>(lineType[0]));
        })

        .def("currentTotal", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.currentTotal(static_cast<unsigned char>(lineType[0]));
        });
}

