#pragma once
#include <vector>

using namespace std;

struct Node {
    double x, y, z;
    int ref_tag; // Tag de référence pour les frontières
};

struct Quad {
    int n1, n2, n3, n4; // Indices des 4 nœuds du quadrilatère
    int ref_tag;
};

class Mesh
{
private:
    double L_x;
    double L_y;
    int n_x;
    int n_y;
public:
    Mesh(double L_x, double L_y, int n_x, int n_y);
    vector<Node> Mesh_points() const;
    vector<Quad> Mesh_quads() const;
    // getters pour les dimensions et le nombre de points
    double get_Lx() const { return L_x; }
    double get_Ly() const { return L_y; }
    int get_nx() const { return n_x; }
    int get_ny() const { return n_y; }
    double get_dx() const { return L_x / (n_x - 1); }
    double get_dy() const { return L_y / (n_y - 1); }
};