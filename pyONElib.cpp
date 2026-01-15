#include <cstring>
#include "ONElib.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
        .def("readLine", [](ONEfile &self) -> bool {
            char result = self.readLine();
            return result != 0;  // Return true if we read a line, false at EOF
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
        // getDNA2bit returns 2-bit compressed DNA as Python bytes
        // 4 bases per byte, so byte_length = (num_bases + 3) / 4
        .def("getDNA2bit", [](ONEfile &self) -> py::bytes {
            uint8_t *data = self.getDNA2bit();
            int64_t num_bases = self.listLength();
            size_t byte_length = (num_bases + 3) / 4;
            return py::bytes(reinterpret_cast<char*>(data), byte_length);
        })
        .def("getString", &ONEfile::getString)
        .def("getStringList", &ONEfile::getStringList)
        .def("getComment", &ONEfile::getComment)
        .def("gotoObject", &ONEfile::gotoObject)
        .def("lineType", &ONEfile::lineType)
        .def("lineNumber", &ONEfile::lineNumber)
        
        // Metadata access methods
        .def("givenCount", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.count(lineType[0]);  // For read files, count() returns given.count
        })

        .def("givenMax", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.max(lineType[0]);  // For read files, max() returns given.max
        })

        .def("givenTotal", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.total(lineType[0]);  // For read files, total() returns given.total
        })

        .def("currentCount", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            // For write files, count() returns accum.count (current)
            return self.count(lineType[0]);
        })

        .def("currentMax", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.max(lineType[0]);
        })

        .def("currentTotal", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.total(lineType[0]);
        })

        // Write methods
        .def("setInt", &ONEfile::setInt)
        .def("setReal", &ONEfile::setReal)
        .def("setChar", &ONEfile::setChar)

        .def("writeLine", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            self.writeLine(lineType[0]);
        })

        .def("writeLine", [](ONEfile &self, const std::string &lineType, const std::string &s) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            self.writeLine(lineType[0], s);
        })

        .def("writeLine", [](ONEfile &self, const std::string &lineType, const std::vector<std::string> &slist) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            self.writeLine(lineType[0], slist);
        })

        .def("writeLineIntList", [](ONEfile &self, const std::string &lineType, const std::vector<int64_t> &data) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            self.writeLine(lineType[0], data.size(), (void*)data.data());
        })

        .def("writeLineRealList", [](ONEfile &self, const std::string &lineType, const std::vector<double> &data) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            self.writeLine(lineType[0], data.size(), (void*)data.data());
        })

        .def("writeLineDNA2bit", [](ONEfile &self, const std::string &lineType, int64_t dnaLen, py::bytes dnaBuf) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            std::string buf = dnaBuf;
            self.writeLineDNA2bit(lineType[0], dnaLen, (uint8_t*)buf.data());
        })

        .def("writeComment", &ONEfile::writeComment)

        // Provenance and reference methods
        .def("addProvenance", &ONEfile::addProvenance)
        .def("inheritProvenance", &ONEfile::inheritProvenance)
        .def("addReference", &ONEfile::addReference)
        .def("inheritReference", &ONEfile::inheritReference)

        // Additional utility methods
        .def("fileName", &ONEfile::fileName)

        .def("count", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.count(lineType[0]);
        })

        .def("max", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.max(lineType[0]);
        })

        .def("total", [](ONEfile &self, const std::string &lineType) {
            if (lineType.length() != 1) {
                throw std::invalid_argument("lineType must be a single character");
            }
            return self.total(lineType[0]);
        })

        .def("checkSchema", &ONEfile::checkSchema);
}

