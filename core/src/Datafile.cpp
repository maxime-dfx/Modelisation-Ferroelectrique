#include "Datafile.h"
#include "Logger.h" 

Datafile::Datafile(const std::string& filename) {
    try {
        m_config = toml::parse_file(filename);
    } catch (const toml::parse_error& err) {
        Logger::error("Erreur lors de la lecture du fichier TOML : " + std::string(err.what()));
        exit(1);
    }
}

