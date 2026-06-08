#include "Datafile.h"
#include <stdexcept>
#include <iostream>

using namespace std;

Datafile::Datafile(const string& filename) {
    try {
        m_config = toml::parse(filename);
    } catch (const std::exception& err) {
        throw runtime_error(string("Echec lecture configuration '") + filename + "' : " + err.what());
    }
}

