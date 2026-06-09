#include "Simulation.h"
#include "Logger.h"
#include "Datafile.h"
#include <iostream>
#include <vector>
#include <fstream>

Simulation::Simulation(PhaseFieldFerro& physics, ResultsExporter& exporter, const Datafile& config)
    : m_physics(physics), m_exporter(exporter), m_current_step(0) {
    
    m_max_steps = config.get_max_steps(); // Récupération du nombre maximum de pas de temps
    m_save_frequency = config.get_save_frequency(); // Récupération de la fréquence d'enregistrement
    m_tolerance = config.get_tolerance(); // Récupération de la tolérance de convergence
    m_output_dir = config.getOutputDir(); // Récupération du répertoire de sortie
    m_vtk_output_dir = config.getVTKOutputDir(); // Récupération du répertoire de sortie pour les fichiers VTK

    m_convergence_data = m_output_dir + "/" + config.getOutputConvergenceDir() + "/" + config.getOutputConvergencePolarization(); // Fichier de données pour la convergence
}

void Simulation::save_current_state() {
    // Fichier unique pour toutes les grandeurs du pas de temps actuel
    std::string filename = m_output_dir + "/" + m_vtk_output_dir + "/state_" + std::to_string(m_current_step) + ".vtk";

    // Les conteneurs qui vont lister ce qu'on exporte
    std::vector<std::string> scalar_names;
    std::vector<const std::vector<double>*> scalar_data;

    std::vector<std::string> vector_names;
    std::vector<const std::vector<double>*> vector_data;

    // Déclaration des tampons mémoires au niveau de la fonction 
    // (Pour garantir que leurs pointeurs restent valides jusqu'à l'appel final)
    std::vector<double> P_combined;
    std::vector<double> E_combined;
    std::vector<double> U_combined;
    std::vector<double> phi_vec;

    // --- 1. Export de la Polarisation (Conditionnel) ---
    if (m_physics.is_thermodynamics_enabled()) {
        const auto& px_data = m_physics.get_Px();
        const auto& py_data = m_physics.get_Py();
        P_combined.reserve(px_data.size() * 2); 
        for(size_t i = 0; i < px_data.size(); ++i) {
            P_combined.push_back(px_data[i]);
            P_combined.push_back(py_data[i]);
        }
        vector_names.push_back("Polarisation");
        vector_data.push_back(&P_combined);
    }

    // --- 2. Export de l'Électrostatique (Conditionnel) ---
    if (m_physics.is_electrostatics_enabled()) {
        // Champ Vectoriel (Ex, Ey)
        const auto& ex_data = m_physics.get_Ex();
        const auto& ey_data = m_physics.get_Ey();
        E_combined.reserve(ex_data.size() * 2);
        for(size_t i = 0; i < ex_data.size(); ++i) {
            E_combined.push_back(ex_data[i]);
            E_combined.push_back(ey_data[i]);
        }
        vector_names.push_back("Champ_Electrique");
        vector_data.push_back(&E_combined);

        // Champ Scalaire (Potentiel)
        const Eigen::VectorXd& phi_data = m_physics.get_potential();
        phi_vec.assign(phi_data.data(), phi_data.data() + phi_data.size());
        scalar_names.push_back("Potentiel");
        scalar_data.push_back(&phi_vec);
    }

    // --- 3. Export de la Mécanique (Conditionnel) ---
    if (m_physics.is_mechanics_enabled()) {
        // Champ Vectoriel (ux, uy)
        const auto& ux_data = m_physics.get_ux();
        const auto& uy_data = m_physics.get_uy();
        U_combined.reserve(ux_data.size() * 2);
        for(size_t i = 0; i < ux_data.size(); ++i) {
            U_combined.push_back(ux_data[i]);
            U_combined.push_back(uy_data[i]);
        }
        vector_names.push_back("Deplacement");
        vector_data.push_back(&U_combined);

        // Champs Scalaires (Contraintes)
        scalar_names.push_back("Sigma_xx");
        scalar_data.push_back(&m_physics.get_sig_xx());
        
        scalar_names.push_back("Sigma_yy");
        scalar_data.push_back(&m_physics.get_sig_yy());
        
        scalar_names.push_back("Sigma_xy");
        scalar_data.push_back(&m_physics.get_sig_xy());
    }

    // --- APPEL UNIQUE D'EXPORTATION ---
    m_exporter.exportMultiPhysicsVTK(filename, scalar_names, scalar_data, vector_names, vector_data);
}

void Simulation::run() {
    bool converged = false;
    m_current_step = 0;
    std::ofstream convergence_log(m_convergence_data);
    convergence_log << "Step\tMaxChange\n"; 
    Logger::info("Debut de la simulation temporelle...");
    save_current_state();

    while (!converged && m_current_step < m_max_steps) {
        m_current_step++;
        double max_change = m_physics.compute_one_step();
        convergence_log << m_current_step << "\t" << max_change << "\n";

        if (max_change < m_tolerance) {
            converged = true;
            Logger::success("Convergence atteinte a l'etape " + std::to_string(m_current_step) + " ! (Variation max = " + std::to_string(max_change) + ")");
        }

        if (m_current_step % m_save_frequency == 0 || converged) {
            save_current_state();
        }   
    }

    if (!converged) {
        Logger::warning("\n[ATTENTION] Simulation arretee par la limite de temps (" + std::to_string(m_max_steps) + " pas) sans converger.");
    }
    convergence_log.close();
}