#include "Validation.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "Logger.h"

// Constante Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void Validation::validation_thermodynamique(const Datafile& datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter) {
    Logger::info("Demarrage de la validation thermodynamique (Norme L2)...");
    
    std::vector<Node> nodes = msh.Mesh_points();
    const std::vector<double>& Py_num = physics.get_Py(); 
    
    // Récupération des paramètres (idéalement liés à ton Datafile)
    double alpha = datafile.get_alpha(); 
    double beta = datafile.get_beta();
    double G = datafile.get_G();
    
    double P0 = std::sqrt(-alpha / (2.0 * beta)); // Polarisation d'équilibre
    double x0 = datafile.get_L()/2; // Position du centre de la paroi
    double epsilon = std::sqrt(-2.0 * G / alpha); // Largeur de la transition
    
    double L2_error = 0.0;
    double L2_norm_exact = 0.0;

    // Vecteurs pour stocker les colonnes de données à exporter
    std::vector<double> vec_x;
    std::vector<double> vec_P_exact;
    std::vector<double> vec_P_num;
    std::vector<double> vec_diff;
    
    vec_x.reserve(nodes.size());
    vec_P_exact.reserve(nodes.size());
    vec_P_num.reserve(nodes.size());
    vec_diff.reserve(nodes.size());

    double y_target = msh.get_ny() * msh.get_dy() / 2.0; 
    double tolerance = msh.get_dy() / 2.0; 

    for (size_t idx = 0; idx < nodes.size(); ++idx) {
        double x = nodes[idx].x;
        double y = nodes[idx].y;
        
        double P_exact = P0 * std::tanh((x - x0) / epsilon); 
        
        double diff = P_exact - Py_num[idx];
        L2_error += diff * diff * msh.get_dx() * msh.get_dy();
        L2_norm_exact += P_exact * P_exact * msh.get_dx() * msh.get_dy();
        
        if (std::abs(y - y_target) <= tolerance) {
            vec_x.push_back(x);
            vec_P_exact.push_back(P_exact);
            vec_P_num.push_back(Py_num[idx]);
            vec_diff.push_back(diff);
        }
    }

    L2_error = std::sqrt(L2_error);
    L2_norm_exact = std::sqrt(L2_norm_exact);
    double relative_error = L2_error / L2_norm_exact;

    Logger::info("Erreur relative L2 (Thermodynamique) : " + std::to_string(relative_error));

    // --- Appel à l'exportation des données ---
    std::vector<std::string> column_names = {"x", "P_exact", "P_num", "Erreur_locale"};
    std::vector<std::vector<double>> data_columns = {vec_x, vec_P_exact, vec_P_num, vec_diff};
    
    std::string export_path = "../data/validation_thermodynamique/thermo_convergence.dat";
    exporter.exportConvergenceData(export_path, column_names, data_columns);
    
    std::vector<double> vec_P_combined;
    vec_P_combined.reserve(vec_P_exact.size() * 2);
    for(size_t i = 0; i < vec_P_exact.size(); ++i) {
        vec_P_combined.push_back(vec_P_exact[i]);
        vec_P_combined.push_back(vec_P_num[i]);
    } 
    exporter.exportStructuredVectorVTK("../data/validation_thermodynamique/VTK/thermo_field.vtk", {"Py"}, vec_P_combined);
    Logger::info("Profil 1D exporte vers : " + export_path);
}

void Validation::validation_electrostatique(const Datafile& datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter) {
    Logger::info("Demarrage de la validation electrostatique (MMS)...");
    
    std::vector<Node> nodes = msh.Mesh_points();
    
    // Paramètres physiques
    double epsilon_0 = datafile.get_epsilon_0(); 
    double epsilon_r = datafile.get_epsilon_r();
    double P0 = 0.26; 
    
    double Lx = msh.get_nx() * msh.get_dx();
    double k = 2.0 * M_PI / Lx; // Nombre d'onde pour une période complète
    
    // 1. Injection de la Solution Manufacturée
    std::vector<double>& Px_mod = const_cast<std::vector<double>&>(physics.get_Px());
    std::vector<double>& Py_mod = const_cast<std::vector<double>&>(physics.get_Py());
    
    for (size_t idx = 0; idx < nodes.size(); ++idx) {
        Px_mod[idx] = P0 * std::sin(k * nodes[idx].x);
        Py_mod[idx] = 0.0; 
    }
    
    // 2. Appel du solveur de Poisson isolément
    physics.solve_electrostatics();
    const Eigen::VectorXd& phi_num = physics.get_potential(); 

    // 3. Comparaison avec la solution analytique
    double L2_error = 0.0;
    double L2_norm_exact = 0.0;
    
    std::vector<double> vec_x;
    std::vector<double> vec_phi_exact;
    std::vector<double> vec_phi_num;
    std::vector<double> vec_diff;

    int j_target = msh.get_ny() / 2; // Capture de la ligne médiane

    for (size_t idx = 0; idx < nodes.size(); ++idx) {
        double x = nodes[idx].x;
        int j = idx / msh.get_nx(); // Coordonnée Y sur la grille
        
        // On aligne la constante d'intégration avec le point d'ancrage (x=0)
        double phi_exact_brut = -(P0 / (k * epsilon_0 * epsilon_r)) * std::cos(k * x);
        double phi_exact_offset = -(P0 / (k * epsilon_0 * epsilon_r)) * std::cos(0.0);
        double phi_exact = phi_exact_brut - phi_exact_offset; 
        
        double diff = phi_exact - phi_num[idx];
        L2_error += diff * diff * msh.get_dx() * msh.get_dy();
        L2_norm_exact += phi_exact * phi_exact * msh.get_dx() * msh.get_dy();
        
        if (j == j_target) {
            vec_x.push_back(x);
            vec_phi_exact.push_back(phi_exact);
            vec_phi_num.push_back(phi_num[idx]);
            vec_diff.push_back(diff);
        }
    }
    
    L2_error = std::sqrt(L2_error);
    L2_norm_exact = std::sqrt(L2_norm_exact);
    double relative_error = L2_error / L2_norm_exact;
    Logger::info("Erreur relative L2 (Electrostatique) : " + std::to_string(relative_error));

    // Exportation des données
    std::vector<std::string> column_names = {"x", "phi_exact", "phi_num", "Erreur_locale"};
    std::vector<std::vector<double>> data_columns = {vec_x, vec_phi_exact, vec_phi_num, vec_diff};
    
    std::string export_path = "../data/validation_electrostatique/electro_convergence.dat";
    exporter.exportConvergenceData(export_path, column_names, data_columns);
    Logger::info("Profil Electrostatique exporte vers : " + export_path);
}

void Validation::validation_fracture(const Datafile& datafile, const Mesh& msh, const PhaseFieldFerro& physics, ResultsExporter& exporter) {
    Logger::info("Demarrage de la validation fracture (Miehe/AT2)...");
    
    /* ATTENTE DE L'IMPLEMENTATION DE LA MECANIQUE DE RUPTURE
    std::vector<Node> nodes = msh.Mesh_points();
    const std::vector<double>& v_num = physics.get_fracture_field(); 
    
    double l = 0.05; 
    double x0 = msh.get_nx() * msh.get_dx() / 2.0; 
    
    double L2_error = 0.0;
    for (size_t idx = 0; idx < nodes.size(); ++idx) {
        double x = nodes[idx].x;
        double v_exact = std::exp(-std::abs(x - x0) / l);
        
        double diff = v_exact - v_num[idx];
        L2_error += diff * diff * msh.get_dx() * msh.get_dy();
    }
    
    L2_error = std::sqrt(L2_error);
    Logger::info("Erreur absolue L2 (Fracture) : " + std::to_string(L2_error));
    */
}

void Validation::run_all_validations(Datafile datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter) {
    if (datafile.enable_validation_thermo() == true) {
        validation_thermodynamique(datafile, msh, physics, exporter);
    }
    if (datafile.enable_validation_electrostatique() == true) {
        validation_electrostatique(datafile, msh, physics, exporter );
    }
    if (datafile.enable_validation_fracture() == true) {
        validation_fracture(datafile, msh, physics, exporter);
    }
}