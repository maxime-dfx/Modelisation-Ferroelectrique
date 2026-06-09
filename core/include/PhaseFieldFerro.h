#pragma once
#include <vector>
#include "Datafile.h"
#include "Mesh.h"
#include "DofConstraint.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>


class PhaseFieldFerro {
private:
    const Mesh& m_mesh;

    double L, alpha, beta, gamma, G, epsilon_0, epsilon_r, dt, dx, dy;
    int nx, ny; 
    std::vector<double> Px, Px_new, Py, Py_new;
    std::vector<double> Ex, Ey;

    Eigen::SparseMatrix<double> A_laplacian; // Matrice de diffusion
    Eigen::SparseMatrix<double> A_imex; // Matrice du schéma IMEX
    Eigen::VectorXd b_charges; // Terme source pour la diffusion
    Eigen::VectorXd  x_potential; // Solution pour la diffusion
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> poisson_solver; // Solveur pour la diffusion
    Eigen::SparseLU<Eigen::SparseMatrix<double>> imex_solver; // Solveur pour le schéma IMEX

    bool m_enable_thermodynamics; // Activer ou désactiver la contribution thermodynamique
    bool m_enable_electrostatics; // Activer ou désactiver la contribution électrostatique
    bool m_enable_mecanics; // Activer ou désactiver la contribution mécanique
    bool m_enable_fracture; // Activer ou désactiver la contribution de la fracture

    bool m_dirichlet_x; // Activer ou désactiver les conditions de Dirichlet sur les bords verticaux
    bool m_dirichlet_y; // Activer ou désactiver les conditions de Dir
    double m_bc_left_value; // Valeur de Dirichlet à gauche
    double m_bc_right_value; // Valeur de Dirichlet à droite
    double m_bc_top_value; // Valeur de Dirichlet en haut
    double m_bc_bottom_value; // Valeur de Dirichlet en bas

    // --- Paramètres Matériaux ---
    double C11, C12, C44;
    double Q11, Q12, Q44;

    // --- Variables Mécaniques ---
    std::vector<double> ux, uy; // Déplacements nodaux
    std::vector<double> sig_xx, sig_yy, sig_xy; // Composantes de la contrainte

    // --- Matrices et Solveurs Mécaniques ---
    Eigen::SparseMatrix<double> K_global; // Matrice de rigidité globale (2N x 2N)
    Eigen::VectorXd F_global;             // Vecteur des forces globales (2N)
    Eigen::VectorXd U_global;             // Solution des déplacements (2N)
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> mechanics_solver; 

    std::vector<DofConstraint> bc_phi; // CL pour le potentiel
    std::vector<DofConstraint> bc_ux; // CL pour u_x
    std::vector<DofConstraint> bc_uy; // CL pour u_y

    // Méthodes pour l'assemblage FEM

    void build_poisson_matrix();
    void build_imex_matrix();
    void compute_electric_field();

    Eigen::Matrix<double, 8, 8> build_Ke();
    Eigen::Matrix<double, 3, 8> get_B_matrix(double xi, double eta);
    void build_mechanics_matrix();
    void solve_mechanics();
    void compute_stresses();

    void solve_fracture();
    double update_polarization();
    void apply_dirichlet_conditions(Eigen::VectorXd& b);

    void init_electrostatics_bcs();
    bool fix_bottom, fix_top, fix_left, fix_right;
    double val_bottom, val_top, val_left, val_right;
    void init_mechanical_bcs();
    bool bottom_x_fixed, bottom_y_fixed;
    bool top_x_fixed, top_y_fixed;
    bool left_x_fixed, left_y_fixed;
    bool right_x_fixed, right_y_fixed;
    double bottom_x_value, bottom_y_value;
    double top_x_value,    top_y_value;
    double left_x_value,   left_y_value;
    double right_x_value,  right_y_value;



public:
    // Signatures uniquement
    PhaseFieldFerro(const Datafile& datafile, const Mesh& mesh);
    void initialize_position();
    void initialize_position_thermodynamique_validation();
    void initialize_position_electrostatique_validation();
    void initialize_position_mecanique_validation();
    void initialize_position_fracture_validation();

    void solve_electrostatics();

    // Les getters d'une ligne restent ici pour la performance (inline)
    const std::vector<double>& get_Px() const { return Px; }
    const std::vector<double>& get_Py() const { return Py; }
    const Eigen::VectorXd& get_potential() const { return x_potential; }

    double compute_one_step();

    bool is_thermodynamics_enabled() const { return m_enable_thermodynamics; }
    bool is_mechanics_enabled() const { return m_enable_mecanics; }

    
    const std::vector<double>& get_ux() const { return ux; }
    const std::vector<double>& get_uy() const { return uy; }
    const std::vector<double>& get_sig_xx() const { return sig_xx; }
    const std::vector<double>& get_sig_yy() const { return sig_yy; }
    const std::vector<double>& get_sig_xy() const { return sig_xy; }

    bool is_electrostatics_enabled() const { return m_enable_electrostatics; }
    const std::vector<double>& get_Ex() const { return Ex; }
    const std::vector<double>& get_Ey() const { return Ey; }
};