#pragma once
#include <string>
#include <iostream>
#include "toml.hpp"

template <typename T>
T get_val(const toml::table& config, const std::string& section, const std::string& key, T default_value) {
    if (auto node = config[section][key]) {
        return node.value<T>().value_or(default_value);
    }
    return default_value;
}

class Datafile {
private:
    toml::table m_config; 

public:
    Datafile(const std::string& filename);

    // --- Paramètres d'Entrée / Sortie (IO) ---
    std::string getMeshCreateFile() const { return get_val<std::string>(m_config, "io", "mesh_create_file", "default.mesh"); }

    std::string getOutputDir() const { return get_val<std::string>(m_config, "io", "output_dir", "data/resultats"); }
    std::string getVTKOutputDir() const { return get_val<std::string>(m_config, "io", "vtk_output_dir", "data/vtk"); }
    std::string getOutputConvergenceDir() const { return get_val<std::string>(m_config, "io", "convergence_output_dir", "convergence"); }

    std::string getOutputPolarizationVtk() const { return get_val<std::string>(m_config, "io", "output_polarization_vtk", "Polarization.vtk"); }
    std::string getOutputConvergencePolarization() const { return get_val<std::string>(m_config, "io", "output_convergence_polarization", "convergence_polarization.dat"); }

    // --- Paramètres de Profiling (Chrono) ---
    bool get_chrono_mesh() const { return get_val<bool>(m_config, "chrono", "temps_calcul_maillage", false); }
    bool get_chrono_export_mesh() const { return get_val<bool>(m_config, "chrono", "temps_export_maillage", false); }
    bool get_chrono_export_VTK_P() const { return get_val<bool>(m_config, "chrono", "temps_export_vtk_p", false); }
    bool get_chrono_global() const { return get_val<bool>(m_config, "chrono", "temps_global", false); }
    bool get_chrono_validation() const { return get_val<bool>(m_config, "chrono", "temps_validation", false); }

    // --- Paramètres du Maillage (Mesh) ---
    double get_L_x() const { return get_val<double>(m_config, "mesh", "L_x", 1.0); }
    double get_L_y() const { return get_val<double>(m_config, "mesh", "L_y", 1.0); }
    int get_n_x() const { return get_val<int>(m_config, "mesh", "n_x", 5); }
    int get_n_y() const { return get_val<int>(m_config, "mesh", "n_y", 5); }

    // --- Paramètres physiques (Physics) ---
    double get_L() const { return get_val<double>(m_config, "physics", "L", 1.0); }
    double get_alpha() const { return get_val<double>(m_config, "physics", "alpha", 1.0); }
    double get_beta() const { return get_val<double>(m_config, "physics", "beta", 1.0); }
    double get_gamma() const { return get_val<double>(m_config, "physics", "gamma", 1.0); }
    double get_G() const { return get_val<double>(m_config, "physics", "G", 1.0); }
    double get_dt() const { return get_val<double>(m_config, "physics", "dt", 0.01); }
    double get_epsilon_0() const { return get_val<double>(m_config, "physics", "epsilon_0", 0.131); }
    double get_epsilon_r() const { return get_val<double>(m_config, "physics", "epsilon_r", 1.0); }

    // --- Paramètres simulation (Simulation) ---
    int get_max_steps() const { return get_val<int>(m_config, "simulation", "max_steps_P", 10000); }
    int get_save_frequency() const { return get_val<int>(m_config, "simulation", "save_frequency_P", 50); }
    double get_tolerance() const { return get_val<double>(m_config, "simulation", "tolerance_P", 1e-6); }
    bool calcul_mesh() const { return get_val<bool>(m_config, "io", "calcul_mesh", true); }

    bool enable_thermodynamics() const { return get_val<bool>(m_config, "physics", "enable_thermodynamics", true); }
    bool enable_electrostatics() const { return get_val<bool>(m_config, "physics", "enable_electrostatics", false); }
    bool enable_mecanics() const { return get_val<bool>(m_config, "physics", "enable_mecanics", false); }
    bool enable_fracture() const { return get_val<bool>(m_config, "physics", "enable_fracture", false); }

    double get_Q11() const { return get_val<double>(m_config, "physics", "Q11", 0.0); }
    double get_Q12() const { return get_val<double>(m_config, "physics", "Q12", 0.0); }
    double get_Q44() const { return get_val<double>(m_config, "physics", "Q44", 0.0); }
    double get_C11() const { return get_val<double>(m_config, "physics", "C11", 0.0); }
    double get_C12() const { return get_val<double>(m_config, "physics", "C12", 0.0); }
    double get_C44() const { return get_val<double>(m_config, "physics", "C44", 0.0); }


    // --- Paramètres de validation (Validation) ---
    bool enable_validation_thermo() const { return get_val<bool>(m_config, "validation", "enable_validation_thermodynamique", false); }
    bool enable_validation_electrostatique() const { return get_val<bool>(m_config, "validation", "enable_validation_electrostatique", false); }
    bool enable_validation_mecanique() const { return get_val<bool>(m_config, "validation", "enable_validation_mecanique", false); }
    bool enable_validation_fracture() const { return get_val<bool>(m_config, "validation", "enable_validation_fracture", false); }

    // --- Paramètres des conditions aux limites Electrostatiques ---
    bool get_fix_bottom_phi() const { return get_val<bool>(m_config, "electrostatics_bcs", "fix_bottom_phi", false); }
    double get_bc_value_bottom_phi() const { return get_val<double>(m_config, "electrostatics_bcs", "value_bottom_phi", 0.0); }
    
    bool get_fix_top_phi() const { return get_val<bool>(m_config, "electrostatics_bcs", "fix_top_phi", false); }
    double get_bc_value_top_phi() const { return get_val<double>(m_config, "electrostatics_bcs", "value_top_phi", 0.0); }
    
    bool get_fix_left_phi() const { return get_val<bool>(m_config, "electrostatics_bcs", "fix_left_phi", false); }
    double get_bc_value_left_phi() const { return get_val<double>(m_config, "electrostatics_bcs", "value_left_phi", 0.0); }
    
    bool get_fix_right_phi() const { return get_val<bool>(m_config, "electrostatics_bcs", "fix_right_phi", false); }
    double get_bc_value_right_phi() const { return get_val<double>(m_config, "electrostatics_bcs", "value_right_phi", 0.0); }

    // --- Paramètres des conditions aux limites Mecaniques ---
    bool enable_dirichlet_x() const { return get_val<bool>(m_config, "boundary_conditions", "enable_dirichlet_x", false); }
    bool enable_dirichlet_y() const { return get_val<bool>(m_config, "boundary_conditions", "enable_dirichlet_y", false); }
    double dirichlet_value_left() const { return get_val<double>(m_config, "boundary_conditions", "dirichlet_left_value", 0.0); }
    double dirichlet_value_right() const { return get_val<double>(m_config, "boundary_conditions", "dirichlet_right_value", 0.0); }
    double dirichlet_value_top() const { return get_val<double>(m_config, "boundary_conditions", "dirichlet_top_value", 0.0); }
    double dirichlet_value_bottom() const { return get_val<double>(m_config, "boundary_conditions", "dirichlet_bottom_value", 0.0); }

    bool get_fix_bottom_x() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_bottom_x", true); }
    double get_bc_value_bottom_x() const { return get_val<double>(m_config, "mechanics_bcs", "value_bottom_x", 0.0); }
    bool get_fix_bottom_y() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_bottom_y", true); }
    double get_bc_value_bottom_y() const { return get_val<double>(m_config, "mechanics_bcs", "value_bottom_y", 0.0); }

    bool get_fix_top_x() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_top_x", false); }
    double get_bc_value_top_x() const { return get_val<double>(m_config, "mechanics_bcs", "value_top_x", 0.0); }
    bool get_fix_top_y() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_top_y", false); }
    double get_bc_value_top_y() const { return get_val<double>(m_config, "mechanics_bcs", "value_top_y", 0.0); }

    bool get_fix_left_x() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_left_x", false); }
    double get_bc_value_left_x() const { return get_val<double>(m_config, "mechanics_bcs", "value_left_x", 0.0); }
    bool get_fix_left_y() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_left_y", false); }
    double get_bc_value_left_y() const { return get_val<double>(m_config, "mechanics_bcs", "value_left_y", 0.0); }

    bool get_fix_right_x() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_right_x", false); }
    double get_bc_value_right_x() const { return get_val<double>(m_config, "mechanics_bcs", "value_right_x", 0.0); }
    bool get_fix_right_y() const { return get_val<bool>(m_config, "mechanics_bcs", "fix_right_y", false); }
    double get_bc_value_right_y() const { return get_val<double>(m_config, "mechanics_bcs", "value_right_y", 0.0); }


};