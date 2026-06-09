#include <iostream>
#include <filesystem>
#include "Mesh.h"
#include "PhaseFieldFerro.h"
#include "ResultsExporter.h"
#include "Simulation.h" 
#include "Chrono.h"
#include "Datafile.h"
#include "Logger.h" 
#include "Validation.h"

int main() {
    try {
        // ==========================================
        // 1. Initialisation et Configuration
        // ==========================================
        Datafile config("../input/config.toml");
        std::filesystem::create_directories(config.getOutputDir());
        std::filesystem::create_directories(config.getVTKOutputDir()); 
        
        Chrono chrono_global, chrono_etape;
        if (config.get_chrono_global()) chrono_global.start();

        Logger::info("Debut du calcul...");

        // ==========================================
        // 2. Calculs et Traitements (Maillage)
        // ==========================================
        if (config.get_chrono_mesh()) chrono_etape.start();

        Mesh msh(config.get_L_x(), config.get_L_y(), 
                 config.get_n_x(), config.get_n_y());

        if (config.get_chrono_mesh()) {
            chrono_etape.stop();
            Logger::time("[MESH Chrono] Generation maillage", chrono_etape.elapsed_ms());
        }

        // ==========================================
        // 3. Entrées / Sorties (Export Maillage)
        // ==========================================
        if (config.get_chrono_export_mesh()) chrono_etape.start();

        ResultsExporter exporter(msh);
        if (!config.calcul_mesh()) {
            exporter.exportToMesh(config.getOutputDir() + "/" + config.getMeshCreateFile() + "_" + std::to_string(config.get_n_x()) + "x" + std::to_string(config.get_n_y()) + ".mesh");
            Logger::warning("Le calcul du maillage est desactive dans le fichier de configuration.");
            exit(1);
        }

        if (config.get_chrono_export_mesh()) {
            chrono_etape.stop();
            Logger::time("[.MESH Chrono] Export disque", chrono_etape.elapsed_ms());
        }

        // ==========================================
        // 4. Physique et Boucle Temporelle
        // ==========================================
        PhaseFieldFerro physics(config, msh);
        Simulation sim(physics, exporter, config);

        // 1. On exécute TOUJOURS la simulation (relaxation du système)
        if (config.get_chrono_export_VTK_P()) chrono_etape.start();

        sim.run();

        if (config.get_chrono_export_VTK_P()) {
            chrono_etape.stop();
            Logger::time("[Polarization Chrono] Simulation & Export", chrono_etape.elapsed_ms());
        }

        if (config.enable_validation_thermo() || config.enable_validation_electrostatique() || config.enable_validation_mecanique() || config.enable_validation_fracture()) {
            if(config.get_chrono_validation()) chrono_etape.start();

            Validation validation(config, msh, physics, exporter);
            validation.run_all_validations();
                
            if(config.get_chrono_validation()) {
                chrono_etape.stop();
                Logger::time("[Validation Chrono] Validations", chrono_etape.elapsed_ms());
            }
        }
        

        // ==========================================
        // 5. Fin et Profiling
        // ==========================================
        if (config.get_chrono_global()) {
            chrono_global.stop();
            Logger::time("TEMPS TOTAL", chrono_global.elapsed_ms());
        } else {
            Logger::success("Execution terminee.");
        }

    } 
    catch (const std::exception& e) {
        Logger::error("Erreur fatale : " + std::string(e.what()));
        return 1; 
    }

    return 0;
}