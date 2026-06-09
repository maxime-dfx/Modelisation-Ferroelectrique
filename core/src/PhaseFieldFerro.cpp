#include "PhaseFieldFerro.h"
#include "Logger.h"
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <algorithm>

// ============================================================
//  Helpers conditions aux limites — polarisation
// ============================================================


static bool is_pol_dirichlet(bool dirichlet_x, bool dirichlet_y,
                              int i, int j, int nx, int ny)
{
    if (dirichlet_x && (i == 0 || i == nx - 1)) return true;
    if (dirichlet_y && (j == 0 || j == ny - 1)) return true;
    return false;
}


static double bc_value_px(bool /*dirichlet_x*/, bool /*dirichlet_y*/,
                           int /*i*/, int /*j*/, int /*nx*/, int /*ny*/,
                           double /*left*/, double /*right*/,
                           double /*bottom*/, double /*top*/)
{
    return 0.0;
}

static double bc_value_py(bool dirichlet_x, bool dirichlet_y,
                           int i, int j, int nx, int ny,
                           double left, double right,
                           double bottom, double top)
{
    // Bords X ont la priorité (électrode principale)
    if (dirichlet_x) {
        if (i == 0)      return left;
        if (i == nx - 1) return right;
    }
    // Bords Y en second (électrode transversale)
    if (dirichlet_y) {
        if (j == 0)      return bottom;
        if (j == ny - 1) return top;
    }
    return 0.0; // fallback (ne devrait pas être atteint)
}

// ============================================================
//  Helpers conditions aux limites — mécanique
// ============================================================

struct EdgeBC {
    bool   x_fixed = false; double x_value = 0.0;
    bool   y_fixed = false; double y_value = 0.0;
};

static void fill_mechanical_bcs(std::vector<DofConstraint>& bc_ux,
                                 std::vector<DofConstraint>& bc_uy,
                                 int nx, int ny,
                                 const EdgeBC& bottom,
                                 const EdgeBC& top,
                                 const EdgeBC& left,
                                 const EdgeBC& right)
{
    // Bord bas (j = 0) → indices 0 … nx-1
    for (int i = 0; i < nx; ++i) {
        if (bottom.x_fixed) { bc_ux[i].is_fixed = true; bc_ux[i].value = bottom.x_value; }
        if (bottom.y_fixed) { bc_uy[i].is_fixed = true; bc_uy[i].value = bottom.y_value; }
    }
    // Bord haut (j = ny-1) → indices (ny-1)*nx … ny*nx-1
    int offset_top = (ny - 1) * nx;
    for (int i = 0; i < nx; ++i) {
        if (top.x_fixed) { bc_ux[offset_top + i].is_fixed = true; bc_ux[offset_top + i].value = top.x_value; }
        if (top.y_fixed) { bc_uy[offset_top + i].is_fixed = true; bc_uy[offset_top + i].value = top.y_value; }
    }
    // Bord gauche (i = 0) → indices j*nx
    for (int j = 0; j < ny; ++j) {
        if (left.x_fixed) { bc_ux[j * nx].is_fixed = true; bc_ux[j * nx].value = left.x_value; }
        if (left.y_fixed) { bc_uy[j * nx].is_fixed = true; bc_uy[j * nx].value = left.y_value; }
    }
    // Bord droit (i = nx-1) → indices j*nx + nx-1
    for (int j = 0; j < ny; ++j) {
        int idx = j * nx + nx - 1;
        if (right.x_fixed) { bc_ux[idx].is_fixed = true; bc_ux[idx].value = right.x_value; }
        if (right.y_fixed) { bc_uy[idx].is_fixed = true; bc_uy[idx].value = right.y_value; }
    }
}

// ============================================================
// Le Constructeur
// ============================================================
PhaseFieldFerro::PhaseFieldFerro(const Datafile& datafile, const Mesh& mesh) : m_mesh(mesh) {
    L = datafile.get_L();
    alpha = datafile.get_alpha();
    beta = datafile.get_beta();
    gamma = datafile.get_gamma();
    G = datafile.get_G();
    dt = datafile.get_dt();
    dx = m_mesh.get_dx();
    dy = m_mesh.get_dy();

    nx = datafile.get_n_x();
    ny = datafile.get_n_y();

    epsilon_0 = datafile.get_epsilon_0();
    epsilon_r = datafile.get_epsilon_r();

    Px.resize(nx * ny, 0.0);
    Py.resize(nx * ny, 0.0);
    Px_new.resize(nx * ny, 0.0);
    Py_new.resize(nx * ny, 0.0);

    int N_nodes = nx * ny;
    Ex.resize(N_nodes, 0.0);
    Ey.resize(N_nodes, 0.0);
    b_charges.resize(N_nodes);
    x_potential.resize(N_nodes);
    A_laplacian.resize(N_nodes, N_nodes);
    A_imex.resize(N_nodes, N_nodes);

    m_enable_thermodynamics = datafile.enable_thermodynamics();
    m_enable_electrostatics = datafile.enable_electrostatics();
    m_enable_mecanics = datafile.enable_mecanics(); 
    m_enable_fracture = datafile.enable_fracture(); 

    m_dirichlet_x = datafile.enable_dirichlet_x();
    m_dirichlet_y = datafile.enable_dirichlet_y();
    if (m_dirichlet_x) {
        m_bc_left_value = datafile.dirichlet_value_left();
        m_bc_right_value = datafile.dirichlet_value_right();
    } else {
        m_bc_left_value = 0.0;
        m_bc_right_value = 0.0;
    }

    if (m_dirichlet_y) {
        m_bc_top_value = datafile.dirichlet_value_top();
        m_bc_bottom_value = datafile.dirichlet_value_bottom();
    } else {
        m_bc_top_value = 0.0;
        m_bc_bottom_value = 0.0;
    }

    // Dans le constructeur de PhaseFieldFerro :
    C11 = datafile.get_C11();
    C12 = datafile.get_C12();
    C44 = datafile.get_C44();
    Q11 = datafile.get_Q11();
    Q12 = datafile.get_Q12();
    Q44 = datafile.get_Q44();

    ux.resize(nx * ny, 0.0);
    uy.resize(nx * ny, 0.0);
    sig_xx.resize(nx * ny, 0.0);
    sig_yy.resize(nx * ny, 0.0);
    sig_xy.resize(nx * ny, 0.0);

    // Taille du système mécanique : 2 * nombre de nœuds
    int N_dof = 2 * nx * ny;
    F_global.resize(N_dof);
    U_global.resize(N_dof);
    K_global.resize(N_dof, N_dof);

    fix_bottom = datafile.get_fix_bottom_phi();
    fix_top    = datafile.get_fix_top_phi();
    fix_left   = datafile.get_fix_left_phi();
    fix_right  = datafile.get_fix_right_phi();   

    val_bottom = fix_bottom ? datafile.get_bc_value_bottom_phi() : 0.0;
    val_top    = fix_top    ? datafile.get_bc_value_top_phi()    : 0.0;
    val_left   = fix_left   ? datafile.get_bc_value_left_phi()   : 0.0;
    val_right  = fix_right  ? datafile.get_bc_value_right_phi()  : 0.0;


    bottom_x_fixed = datafile.get_fix_bottom_x();
    bottom_y_fixed = datafile.get_fix_bottom_y();
    top_x_fixed    = datafile.get_fix_top_x();
    top_y_fixed    = datafile.get_fix_top_y();
    left_x_fixed   = datafile.get_fix_left_x();
    left_y_fixed   = datafile.get_fix_left_y();
    right_x_fixed  = datafile.get_fix_right_x();
    right_y_fixed  = datafile.get_fix_right_y();

    // Valeurs imposées (déplacement ou force) pour chaque bord
    bottom_x_value = bottom_x_fixed ? datafile.get_bc_value_bottom_x() : 0.0;
    bottom_y_value = bottom_y_fixed ? datafile.get_bc_value_bottom_y() : 0.0;
    top_x_value    = top_x_fixed    ? datafile.get_bc_value_top_x()    : 0.0;
    top_y_value    = top_y_fixed    ? datafile.get_bc_value_top_y()    : 0.0;
    left_x_value   = left_x_fixed   ? datafile.get_bc_value_left_x()   : 0.0;
    left_y_value   = left_y_fixed   ? datafile.get_bc_value_left_y()   : 0.0;
    right_x_value  = right_x_fixed  ? datafile.get_bc_value_right_x()  : 0.0;
    right_y_value  = right_y_fixed  ? datafile.get_bc_value_right_y()  : 0.0;



    // ==========================================
    // Initialisation conditionnelle de la physique
    // ==========================================
    if (m_enable_thermodynamics) {
        if (datafile.enable_validation_thermo() == true) {
            initialize_position_thermodynamique_validation();
        } else {
            initialize_position(); 
        }
    } else {
        std::fill(Px.begin(), Px.end(), 0.0);
        std::fill(Py.begin(), Py.end(), 0.0);
    }

    Logger::info("Initialisation de la physique termodynamique : " + std::string(m_enable_thermodynamics ? "ON" : "OFF"));
    Logger::info("Initialisation de la physique electrostatique : " + std::string(m_enable_electrostatics ? "ON" : "OFF"));
    Logger::info("Initialisation de la physique mecanique : " + std::string(m_enable_mecanics ? "ON" : "OFF"));
    Logger::info("Initialisation de la physique fracture : " + std::string(m_enable_fracture ? "ON" : "OFF"));

    init_electrostatics_bcs();
    init_mechanical_bcs();
    if (m_enable_mecanics) build_mechanics_matrix();  
    if (m_enable_electrostatics) build_poisson_matrix();
    if (m_enable_thermodynamics) build_imex_matrix();
}

void PhaseFieldFerro::initialize_position() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    for (size_t i = 0; i < Px.size(); ++i) {
        Px[i] = 0.1 * ((std::rand() / (double)RAND_MAX) - 0.5); 
        Py[i] = 0.1 * ((std::rand() / (double)RAND_MAX) - 0.5);
    }
}

void PhaseFieldFerro::build_poisson_matrix() {
    int N_nodes = nx * ny;   
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(N_nodes * 5);

    double inv_dx2 = 1.0 / (dx * dx);
    double inv_dy2 = 1.0 / (dy * dy);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;
            double diag = 0.0; 

            if (i > 0)      { triplets.push_back(Eigen::Triplet<double>(idx, idx - 1, -inv_dx2)); diag += inv_dx2; }
            if (i < nx - 1) { triplets.push_back(Eigen::Triplet<double>(idx, idx + 1, -inv_dx2)); diag += inv_dx2; }
            if (j > 0)      { triplets.push_back(Eigen::Triplet<double>(idx, idx - nx, -inv_dy2)); diag += inv_dy2; }
            if (j < ny - 1) { triplets.push_back(Eigen::Triplet<double>(idx, idx + nx, -inv_dy2)); diag += inv_dy2; }

            for (int i = 0; i < nx * ny; ++i) {
                if (bc_phi[i].is_fixed) {
                    triplets.push_back(Eigen::Triplet<double>(i, i, 1e30));
                }
            }

            triplets.push_back(Eigen::Triplet<double>(idx, idx, diag));
        }
    }
    
    A_laplacian.setFromTriplets(triplets.begin(), triplets.end());
    poisson_solver.compute(A_laplacian);
    if (poisson_solver.info() != Eigen::Success) { Logger::error("Echec Poisson"); exit(1); }
}

void PhaseFieldFerro::build_imex_matrix() { 
    int N_nodes = nx * ny;   
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(N_nodes * 5);

    double inv_dx2 = 1.0 / (dx * dx);
    double inv_dy2 = 1.0 / (dy * dy);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;

            if (is_pol_dirichlet(m_dirichlet_x, m_dirichlet_y, i, j, nx, ny)) {
                // Imposition forte : équation identité 1·P_new = BC_value
                triplets.push_back(Eigen::Triplet<double>(idx, idx, 1.0));
            } 
            else {
                // Diffusion implicite standard (Ginzburg-Landau / Neumann naturel)
                double diag = 1.0; 
                if (i > 0)      { triplets.push_back({idx, idx - 1,  -dt*L*G*inv_dx2}); diag += dt*L*G*inv_dx2; }
                if (i < nx - 1) { triplets.push_back({idx, idx + 1,  -dt*L*G*inv_dx2}); diag += dt*L*G*inv_dx2; }
                if (j > 0)      { triplets.push_back({idx, idx - nx, -dt*L*G*inv_dy2}); diag += dt*L*G*inv_dy2; }
                if (j < ny - 1) { triplets.push_back({idx, idx + nx, -dt*L*G*inv_dy2}); diag += dt*L*G*inv_dy2; }
                triplets.push_back({idx, idx, diag});
            }
        }
    }
    A_imex.setFromTriplets(triplets.begin(), triplets.end());
    imex_solver.compute(A_imex);
    if (imex_solver.info() != Eigen::Success) { Logger::error("Echec IMEX"); exit(1); }
}

void PhaseFieldFerro::solve_electrostatics() {

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;
            
            // 1. Calcul de la divergence par différences finies centrales
            double px_right = (i < nx - 1) ? Px[idx + 1] : Px[idx];
            double px_left  = (i > 0)      ? Px[idx - 1] : Px[idx];
            double py_up    = (j < ny - 1) ? Py[idx + nx] : Py[idx];
            double py_down  = (j > 0)      ? Py[idx - nx] : Py[idx];
            
            double div_P = (px_right - px_left) / (2.0 * dx) + (py_up - py_down) / (2.0 * dy);
            double total_charge = -div_P;

            // 2. Charges de surface sur les conditions aux limites (murs)
            if (i == 0)      total_charge -= Px[idx] / dx; 
            if (i == nx - 1) total_charge += Px[idx] / dx; 
            if (j == 0)      total_charge -= Py[idx] / dy; 
            if (j == ny - 1) total_charge += Py[idx] / dy; 

            // 3. Remplissage du second membre pour le solveur de Poisson
            b_charges[idx] = total_charge / (epsilon_0 * epsilon_r); 
            
            for (int i = 0; i < nx * ny; ++i) {
                if (bc_phi[i].is_fixed) {
                    b_charges[i] = bc_phi[i].value * 1e30;
                }
            }
        }
    }
    
    x_potential = poisson_solver.solve(b_charges);
    if (poisson_solver.info() != Eigen::Success) { 
        Logger::error("Echec resolution Poisson"); 
        exit(1); 
    }
}

void PhaseFieldFerro::init_electrostatics_bcs() {
    bc_phi.assign(nx * ny, {false, 0.0});

    // Bord bas et haut
    for (int i = 0; i < nx; ++i) {
        if (fix_bottom) { bc_phi[i].is_fixed = true; bc_phi[i].value = val_bottom; }
        if (fix_top)    { bc_phi[(ny - 1) * nx + i].is_fixed = true; bc_phi[(ny - 1) * nx + i].value = val_top; }
    }
    // Bord gauche et droit
    for (int j = 0; j < ny; ++j) {
        if (fix_left)  { bc_phi[j * nx].is_fixed = true; bc_phi[j * nx].value = val_left; }
        if (fix_right) { bc_phi[j * nx + nx - 1].is_fixed = true; bc_phi[j * nx + nx - 1].value = val_right; }
    }

    // Sécurité mathématique absolue : si l'utilisateur ne fixe AUCUN potentiel,
    bool has_dirichlet = false;
    for (const auto& bc : bc_phi) { if (bc.is_fixed) has_dirichlet = true; }
    if (!has_dirichlet) {
        bc_phi[0].is_fixed = true;
        bc_phi[0].value = 0.0;
    }
}

void PhaseFieldFerro::init_mechanical_bcs() {
    bc_ux.assign(nx * ny, {false, 0.0});
    bc_uy.assign(nx * ny, {false, 0.0});

    EdgeBC bottom { bottom_x_fixed, bottom_x_value, bottom_y_fixed, bottom_y_value };
    EdgeBC top    { top_x_fixed,    top_x_value,    top_y_fixed,    top_y_value    };
    EdgeBC left   { left_x_fixed,   left_x_value,   left_y_fixed,   left_y_value   };
    EdgeBC right  { right_x_fixed,  right_x_value,  right_y_fixed,  right_y_value  };

    fill_mechanical_bcs(bc_ux, bc_uy, nx, ny, bottom, top, left, right);
}

Eigen::Matrix<double, 3, 8> PhaseFieldFerro::get_B_matrix(double xi, double eta) {
    Eigen::Matrix<double, 3, 8> B = Eigen::Matrix<double, 3, 8>::Zero();
    
    // Dérivées par rapport à xi et eta
    double dN1_dxi = -0.25 * (1.0 - eta); double dN1_deta = -0.25 * (1.0 - xi);
    double dN2_dxi =  0.25 * (1.0 - eta); double dN2_deta = -0.25 * (1.0 + xi);
    double dN3_dxi =  0.25 * (1.0 + eta); double dN3_deta =  0.25 * (1.0 + xi);
    double dN4_dxi = -0.25 * (1.0 + eta); double dN4_deta =  0.25 * (1.0 - xi);

    double invJ11 = 2.0 / dx; double invJ22 = 2.0 / dy;

    // Remplissage optimisé
    B(0, 0) = dN1_dxi * invJ11; B(1, 1) = dN1_deta * invJ22; B(2, 0) = B(1, 1); B(2, 1) = B(0, 0);
    B(0, 2) = dN2_dxi * invJ11; B(1, 3) = dN2_deta * invJ22; B(2, 2) = B(1, 3); B(2, 3) = B(0, 2);
    B(0, 4) = dN3_dxi * invJ11; B(1, 5) = dN3_deta * invJ22; B(2, 4) = B(1, 5); B(2, 5) = B(0, 4);
    B(0, 6) = dN4_dxi * invJ11; B(1, 7) = dN4_deta * invJ22; B(2, 6) = B(1, 7); B(2, 7) = B(0, 6);

    return B;
}

Eigen::Matrix<double, 8, 8> PhaseFieldFerro::build_Ke() {
    Eigen::Matrix<double, 8, 8> Ke = Eigen::Matrix<double, 8, 8>::Zero();
    Eigen::Matrix<double, 3, 3> D = Eigen::Matrix<double, 3, 3>::Zero();
    D(0, 0) = C11; D(0, 1) = C12; D(1, 0) = C12; D(1, 1) = C11; D(2, 2) = C44;

    double det_J = (dx / 2.0) * (dy / 2.0);
    double pt = 1.0 / std::sqrt(3.0);
    std::vector<std::pair<double, double>> gauss_points = {{-pt, -pt}, {pt, -pt}, {pt, pt}, {-pt, pt}};

    for (const auto& point : gauss_points) {
        Eigen::Matrix<double, 3, 8> B = get_B_matrix(point.first, point.second);
        Ke += B.transpose() * D * B * det_J;
    }
    return Ke;
}

void PhaseFieldFerro::build_mechanics_matrix() {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nx * ny * 64); 

    Eigen::Matrix<double, 8, 8> Ke = build_Ke(); 
    std::vector<Quad> quads = m_mesh.Mesh_quads();

    for (size_t e = 0; e < quads.size(); ++e) {
        const Quad& quad = quads[e];

        int n1 = quad.n1 - 1;
        int n2 = quad.n2 - 1;
        int n3 = quad.n3 - 1;
        int n4 = quad.n4 - 1;

        std::vector<int> dof_indices = {
            2 * n1,     2 * n1 + 1,
            2 * n2,     2 * n2 + 1,
            2 * n3,     2 * n3 + 1,
            2 * n4,     2 * n4 + 1
        };

        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                triplets.push_back(Eigen::Triplet<double>(dof_indices[i], dof_indices[j], Ke(i, j)));
            }
        }
    }

    for (int i = 0; i < nx * ny; ++i) {
        if (bc_ux[i].is_fixed) triplets.push_back(Eigen::Triplet<double>(2*i,   2*i,   1e30));
        if (bc_uy[i].is_fixed) triplets.push_back(Eigen::Triplet<double>(2*i+1, 2*i+1, 1e30));
    }

    K_global.setFromTriplets(triplets.begin(), triplets.end());
    mechanics_solver.compute(K_global);
    if (mechanics_solver.info() != Eigen::Success) { Logger::error("Echec factorisation mecanique"); exit(1); }
}

void PhaseFieldFerro::solve_mechanics() {
    // 1. Remise à zéro du vecteur force à chaque pas de temps
    F_global.setZero(); 

    // 2. Initialisation de la matrice de rigidité D
    Eigen::Matrix<double, 3, 3> D = Eigen::Matrix<double, 3, 3>::Zero();
    D(0, 0) = C11; D(0, 1) = C12; 
    D(1, 0) = C12; D(1, 1) = C11; 
    D(2, 2) = C44;

    // Constantes d'intégration de Gauss (4 points)
    const double det_J = (dx / 2.0) * (dy / 2.0);
    const double pt = 1.0 / std::sqrt(3.0);
    const std::vector<std::pair<double, double>> gauss_points = {
        {-pt, -pt}, {pt, -pt}, {pt, pt}, {-pt, pt}
    };

    const std::vector<Quad>& quads = m_mesh.Mesh_quads();

    // 3. Boucle sur les éléments : Calcul des forces d'électrostriction
    for (const auto& quad : quads) {
        int n1 = quad.n1 - 1; 
        int n2 = quad.n2 - 1;
        int n3 = quad.n3 - 1; 
        int n4 = quad.n4 - 1;

        // Utilisation de std::array (plus rapide que std::vector car allocation statique en mémoire)
        std::array<int, 8> dof_indices = {
            2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n3, 2*n3+1, 2*n4, 2*n4+1
        };

        // Extraction de la polarisation locale sur les nœuds
        double px1 = Px[n1], py1 = Py[n1];
        double px2 = Px[n2], py2 = Py[n2];
        double px3 = Px[n3], py3 = Py[n3];
        double px4 = Px[n4], py4 = Py[n4];

        Eigen::VectorXd Fe = Eigen::VectorXd::Zero(8);

        for (const auto& point : gauss_points) {
            double xi = point.first; 
            double eta = point.second;

            // Fonctions de forme pour interpoler P au point de Gauss
            double N1 = 0.25 * (1.0 - xi) * (1.0 - eta);
            double N2 = 0.25 * (1.0 + xi) * (1.0 - eta);
            double N3 = 0.25 * (1.0 + xi) * (1.0 + eta);
            double N4 = 0.25 * (1.0 - xi) * (1.0 + eta);

            double P_x_gauss = N1*px1 + N2*px2 + N3*px3 + N4*px4;
            double P_y_gauss = N1*py1 + N2*py2 + N3*py3 + N4*py4;

            // Vecteur de déformation propre (e0)
            Eigen::Vector3d e_0;
            e_0(0) = Q11 * P_x_gauss * P_x_gauss + Q12 * P_y_gauss * P_y_gauss;
            e_0(1) = Q12 * P_x_gauss * P_x_gauss + Q11 * P_y_gauss * P_y_gauss;
            e_0(2) = Q44 * P_x_gauss * P_y_gauss;

            Eigen::Matrix<double, 3, 8> B = get_B_matrix(xi, eta);
            
            // OPTIMISATION EIGEN : 
            // 1. Les parenthèses (D * e_0) forcent Eigen à multiplier un 3x3 avec un 3x1 d'abord (très rapide), 
            //    puis de multiplier B^T avec le résultat.
            // 2. noalias() empêche Eigen de créer une matrice temporaire inutile.
            Fe.noalias() += B.transpose() * (D * e_0) * det_J;
        }

        // Assemblage des forces dans le vecteur global
        for (int i = 0; i < 8; ++i) {
            F_global(dof_indices[i]) += Fe(i);
        }
    }

    // 4. Gestion de la méthode de pénalité pour les conditions de Dirichlet
    for (int i = 0; i < nx * ny; ++i) {
        if (bc_ux[i].is_fixed) F_global(2*i)   = bc_ux[i].value * 1e30;
        if (bc_uy[i].is_fixed) F_global(2*i+1) = bc_uy[i].value * 1e30;
    }

    // 5. Résolution du système linéaire
    U_global = mechanics_solver.solve(F_global);

    // 6. Mise à jour des déplacements pour l'exportation
    for (int i = 0; i < nx * ny; ++i) {
        ux[i] = U_global(2 * i);
        uy[i] = U_global(2 * i + 1);
    }
}

void PhaseFieldFerro::compute_stresses() {
    // 1. Initialisation à zéro et création d'un compteur d'adjacence
    std::fill(sig_xx.begin(), sig_xx.end(), 0.0);
    std::fill(sig_yy.begin(), sig_yy.end(), 0.0);
    std::fill(sig_xy.begin(), sig_xy.end(), 0.0);
    std::vector<int> node_counts(nx * ny, 0);

    // Matrice d'élasticité
    Eigen::Matrix<double, 3, 3> D = Eigen::Matrix<double, 3, 3>::Zero();
    D(0, 0) = C11; D(0, 1) = C12; D(1, 0) = C12; D(1, 1) = C11; D(2, 2) = C44;

    // La matrice B évaluée exactement au CENTRE de l'élément (xi = 0, eta = 0)
    Eigen::Matrix<double, 3, 8> B_center = get_B_matrix(0.0, 0.0);
    
    std::vector<Quad> quads = m_mesh.Mesh_quads();

    // 2. Boucle sur les ÉLÉMENTS (et non plus sur les nœuds)
    for (const auto& quad : quads) {
        int n1 = quad.n1 - 1;
        int n2 = quad.n2 - 1;
        int n3 = quad.n3 - 1;
        int n4 = quad.n4 - 1;

        // Extraction des déplacements de cet élément
        Eigen::VectorXd Ue(8);
        Ue << ux[n1], uy[n1], ux[n2], uy[n2],
              ux[n3], uy[n3], ux[n4], uy[n4];

        // Calcul exact de la déformation totale au centre
        Eigen::Vector3d eps_tot = B_center * Ue;

        // Interpolation de la polarisation au centre de l'élément
        double px_c = 0.25 * (Px[n1] + Px[n2] + Px[n3] + Px[n4]);
        double py_c = 0.25 * (Py[n1] + Py[n2] + Py[n3] + Py[n4]);

        // Déformation propre (Électrostriction)
        Eigen::Vector3d eps_0;
        eps_0(0) = Q11 * px_c * px_c + Q12 * py_c * py_c;
        eps_0(1) = Q12 * px_c * px_c + Q11 * py_c * py_c;
        eps_0(2) = Q44 * px_c * py_c;

        // Contrainte constante au centre de l'élément via la Loi de Hooke
        Eigen::Vector3d sig_elem = D * (eps_tot - eps_0);

        // 3. Distribution de la contrainte aux 4 nœuds
        std::vector<int> nodes = {n1, n2, n3, n4};
        for (int n : nodes) {
            sig_xx[n] += sig_elem(0);
            sig_yy[n] += sig_elem(1);
            sig_xy[n] += sig_elem(2);
            node_counts[n]++;
        }
    }

    // 4. Lissage final (On moyenne par le nombre d'éléments qui touchent chaque nœud)
    for (int i = 0; i < nx * ny; ++i) {
        if (node_counts[i] > 0) {
            sig_xx[i] /= node_counts[i];
            sig_yy[i] /= node_counts[i];
            sig_xy[i] /= node_counts[i];
        }
    }
}

void PhaseFieldFerro::solve_fracture() {
    // Placeholder pour la contribution de fracture
    // On pourrait introduire une variable de champ de fracture phi qui modifie localement les coefficients alpha, beta, gamma, ou même introduire une dynamique de champ de fracture qui évolue en parallèle avec la polarisation.
}

void PhaseFieldFerro::compute_electric_field() {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;
            
            if (i == 0)           Ex[idx] = -(x_potential[idx + 1] - x_potential[idx]) / dx;
            else if (i == nx - 1) Ex[idx] = -(x_potential[idx] - x_potential[idx - 1]) / dx;
            else                  Ex[idx] = -(x_potential[idx + 1] - x_potential[idx - 1]) / (2.0 * dx);

            if (j == 0)           Ey[idx] = -(x_potential[idx + nx] - x_potential[idx]) / dy;
            else if (j == ny - 1) Ey[idx] = -(x_potential[idx] - x_potential[idx - nx]) / dy;
            else                  Ey[idx] = -(x_potential[idx + nx] - x_potential[idx - nx]) / (2.0 * dy);
        }
    }
}

double PhaseFieldFerro::compute_one_step() {
    double max_delta = 0.0; 

    if (m_enable_electrostatics) {
        solve_electrostatics();
        compute_electric_field();
    } else {
        std::fill(Ex.begin(), Ex.end(), 0.0);
        std::fill(Ey.begin(), Ey.end(), 0.0);
    }

    if (m_enable_mecanics) {
        solve_mechanics();
        compute_stresses(); 
    }

    if (m_enable_fracture) {
        solve_fracture();
    }

    if (m_enable_thermodynamics) {
        max_delta = update_polarization();
    } else {
        max_delta = 0.0; 
    }

    return max_delta;
}

double PhaseFieldFerro::update_polarization() {
    Eigen::VectorXd b_imex_x(nx * ny);
    Eigen::VectorXd b_imex_y(nx * ny);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j * nx + i;

            if (is_pol_dirichlet(m_dirichlet_x, m_dirichlet_y, i, j, nx, ny)) {
                b_imex_x[idx] = bc_value_px(m_dirichlet_x, m_dirichlet_y, i, j, nx, ny,
                                             0.0, 0.0, 0.0, 0.0);
                b_imex_y[idx] = bc_value_py(m_dirichlet_x, m_dirichlet_y, i, j, nx, ny,
                                             m_bc_left_value, m_bc_right_value,
                                             m_bc_bottom_value, m_bc_top_value);
                continue;
            }

            const double px = Px[idx];
            const double py = Py[idx];

            // --- Landau (isotrope + anisotrope) ---
            const double P2    = px*px + py*py;
            const double iso_x = alpha*px + beta*P2*px;
            const double iso_y = alpha*py + beta*P2*py;
            const double ani_x = gamma*py*py*px;
            const double ani_y = gamma*px*px*py;

            // --- Couplage électrostrictif ---
            double mec_x = 0.0, mec_y = 0.0;
            if (m_enable_mecanics) {
                mec_x = 2.0*(Q11*sig_xx[idx] + Q12*sig_yy[idx])*px
                       + Q44*sig_xy[idx]*py;
                mec_y = 2.0*(Q12*sig_xx[idx] + Q11*sig_yy[idx])*py
                       + Q44*sig_xy[idx]*px;
            }

            // --- Second membre TDGL complet ---
            b_imex_x[idx] = px - dt*L*(iso_x + ani_x - Ex[idx] - mec_x);
            b_imex_y[idx] = py - dt*L*(iso_y + ani_y - Ey[idx] - mec_y);
        }
    }

    Eigen::VectorXd Px_sol = imex_solver.solve(b_imex_x);
    Eigen::VectorXd Py_sol = imex_solver.solve(b_imex_y);

    double max_delta = 0.0;
    for (int i = 0; i < nx * ny; ++i) {
        Px_new[i] = Px_sol[i];
        Py_new[i] = Py_sol[i];
        max_delta = std::max(max_delta, std::abs(Px_new[i] - Px[i]));
        max_delta = std::max(max_delta, std::abs(Py_new[i] - Py[i]));
    }

    std::swap(Px, Px_new);
    std::swap(Py, Py_new);

    return max_delta;
}

void PhaseFieldFerro::initialize_position_thermodynamique_validation() {
    double P0 = std::sqrt(-alpha / (2.0 * beta)); 
    std::vector<Node> nodes = m_mesh.Mesh_points();
    for (size_t idx = 0; idx < nodes.size(); ++idx) {
        Px[idx] = 0.0; 
        int i = idx % nx; 
        if (i < nx / 2) {
            Py[idx] = -P0;
        } else {
            Py[idx] = P0;
        }
    }
}