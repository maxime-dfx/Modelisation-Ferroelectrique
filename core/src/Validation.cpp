#include "Validation.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "Logger.h"

// Constante Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Validation::Validation(const Datafile& datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter)
    : datafile(datafile), msh(msh), physics(physics), exporter(exporter) {}

void Validation::validation_thermodynamique() {
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

void Validation::validation_electrostatique() {
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

void Validation::validation_mecanique() {
    Logger::info("Demarrage de la validation mecanique (Traction uniaxiale)...");
    
    std::vector<Node> nodes = msh.Mesh_points();
    const std::vector<double>& ux_num = physics.get_ux();
    const std::vector<double>& uy_num = physics.get_uy();

    // 1. Paramètres physiques et géométriques
    double Ly = msh.get_Ly();
    double u_top = datafile.get_bc_value_top_y();
    double eps_yy = u_top / Ly;

    double C11 = datafile.get_C11();
    double C12 = datafile.get_C12();

    // Effet Poisson pour le déplacement en X
    double nu_eff = (C11 != 0.0) ? (C12 / C11) : 0.0;
    double eps_xx = -nu_eff * eps_yy;

    // 2. Calcul des erreurs
    double L2_error_u = 0.0;
    double L2_norm_exact = 0.0;

    std::vector<double> vec_y;
    std::vector<double> vec_uy_exact;
    std::vector<double> vec_uy_num;
    std::vector<double> vec_diff;

    size_t nx = static_cast<size_t>(msh.get_nx());
    size_t col_center = nx / 2; 

    for (size_t idx = 0; idx < nodes.size(); ++idx) {
        double x = nodes[idx].x;
        double y = nodes[idx].y;

        // Solution analytique
        double uy_exact = eps_yy * y;
        double ux_exact = eps_xx * x;

        double diff_x = ux_exact - ux_num[idx];
        double diff_y = uy_exact - uy_num[idx];

        // Intégration L2
        L2_error_u += (diff_x * diff_x + diff_y * diff_y) * msh.get_dx() * msh.get_dy();
        L2_norm_exact += (ux_exact * ux_exact + uy_exact * uy_exact) * msh.get_dx() * msh.get_dy();

        // Enregistrement du profil vertical
        if (idx % nx == col_center) {
            vec_y.push_back(y);
            vec_uy_exact.push_back(uy_exact);
            vec_uy_num.push_back(uy_num[idx]);
            vec_diff.push_back(diff_y);
        }
    }

    double relative_error = (L2_norm_exact > 1e-14) ? std::sqrt(L2_error_u) / std::sqrt(L2_norm_exact) : std::sqrt(L2_error_u);
    Logger::info("Erreur relative L2 globale (Deplacements) : " + std::to_string(relative_error));

    // 3. Exportation des données pour tracé Python/Gnuplot
    std::vector<std::string> column_names = {"y", "uy_exact", "uy_num", "Erreur_locale"};
    std::vector<std::vector<double>> data_columns = {vec_y, vec_uy_exact, vec_uy_num, vec_diff};
    
    std::string export_path = "../data/validation_mecanique/meca_convergence.dat";
    exporter.exportConvergenceData(export_path, column_names, data_columns);
    Logger::info("Profil Mecanique 1D exporte vers : " + export_path);
}

void Validation::validation_fracture() {
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

void Validation::run_all_validations() {
    Logger::info("Lancement de toutes les validations...");
    Logger::info("Validation Thermodynamique : " + std::string(datafile.enable_validation_thermo() ? "ON" : "OFF"));
    Logger::info("Validation Electrostatique : " + std::string(datafile.enable_validation_electrostatique() ? "ON" : "OFF"));
    Logger::info("Validation Mecanique : " + std::string(datafile.enable_validation_mecanique() ? "ON" : "OFF"));
    Logger::info("Validation Fracture : " + std::string(datafile.enable_validation_fracture() ? "ON" : "OFF"));
    if (datafile.enable_validation_thermo() == true) {
        validation_thermodynamique();
    }
    if (datafile.enable_validation_electrostatique() == true) {
        validation_electrostatique();
    }
    if (datafile.enable_validation_mecanique() == true) {
        validation_mecanique();
    }
    if (datafile.enable_validation_fracture() == true) {
        validation_fracture();
    }
}