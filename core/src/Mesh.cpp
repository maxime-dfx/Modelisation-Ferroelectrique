#include "Mesh.h"
#include <vector>

using namespace std;

Mesh::Mesh(double L_x, double L_y, int n_x, int n_y) : L_x(L_x), L_y(L_y), n_x(n_x), n_y(n_y) {}

vector<Node> Mesh::Mesh_points() const
{
    vector<Node> mesh_points;
    double dx = get_dx();
    double dy = get_dy();

    mesh_points.reserve(get_nx() * get_ny()); // Réserve de la mémoire pour les nœuds 
    for (int j = 0; j < get_ny(); ++j) {
        for (int i = 0; i < get_nx(); ++i) {
            double x = i * dx;
            double y = j * dy;
            mesh_points.push_back({x, y, 0.0, 0}); // z = 0.0 for a 2D mesh, ref_tag = 0
        }
    }
    return mesh_points;
};

vector<Quad> Mesh::Mesh_quads() const   
{
    vector<Quad> quads;
    quads.reserve((get_nx() - 1) * (get_ny() - 1)); // Réserve de la mémoire pour les quadrilatères

    for (int j = 0; j < get_ny() - 1; ++j) {
        for (int i = 0; i < get_nx() - 1; ++i) {
            int n1 = (j * get_nx()) + i + 1; // Indice du nœud en bas à gauche
            int n2 = n1 + 1;             // Indice du nœud en bas à droite
            int n3 = ((j + 1) * get_nx()) + i + 2; // Indice du nœud en haut à droite
            int n4 = ((j + 1) * get_nx()) + i + 1; // Indice du nœud en haut à gauche
            quads.push_back({n1, n2, n3, n4, 0}); // ref_tag = 0 pour l'intérieur du domaine
        }
    }
    return quads;
};