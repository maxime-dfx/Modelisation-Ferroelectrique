#include <fstream>
#include <iostream>
#include <filesystem>
#include "ResultsExporter.h"
#include "Logger.h" 


ResultsExporter::ResultsExporter(const Mesh& msh) : m_mesh(msh) {}


void ResultsExporter::ensure_directory_exists(const std::string& path) {
    std::filesystem::path filepath(path);
    std::filesystem::path dir = filepath.parent_path();

    if (!dir.empty() && !std::filesystem::exists(dir)) {
        if (!std::filesystem::create_directories(dir)) {
            Logger::error("Echec de la creation du repertoire : " + dir.string());
            exit(1);
        }
    }
}

void ResultsExporter::exportToMesh(const std::string& filename) {
    Logger::info("Demarrage de l'export du maillage vers : " + filename);

    // Vérification et création du dossier
    ensure_directory_exists(filename);

    // 1. Génération des données
    std::vector<Node> nodes = m_mesh.Mesh_points();
    std::vector<Quad> elements = m_mesh.Mesh_quads();

    // 2. Ouverture du fichier
    std::ofstream out(filename); 
    
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    // 3. Écriture des données au format .mesh
    out << "MeshVersionFormatted 1\nDimension 3\n\n";

    out << "Vertices\n" << nodes.size() << "\n";
    for (const auto& node : nodes) {
        out << node.x << " " << node.y << " " << node.z << " " << node.ref_tag << "\n";
    }
    out << "\n";

    out << "Quadrilaterals\n" << elements.size() << "\n";
    for (const auto& quad : elements) {
        out << quad.n1 << " " << quad.n2 << " " << quad.n3 << " " << quad.n4 << " " << quad.ref_tag << "\n";
    }
    out << "\nEnd\n";
    
    out.close();
}

void ResultsExporter::exportStructuredScalarVTK(const std::string& filename, 
                                                const std::vector<std::string>& field_names, 
                                                const std::vector<const std::vector<double>*>& data_fields) {


    // Vérification et création du dossier
    ensure_directory_exists(filename);

    std::ofstream out(filename);
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    // En-tête VTK
    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro Simulation\nASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    
    // Ici on définit la géométrie : nx points en largeur, ny en hauteur, 1 en profondeur
    out << "DIMENSIONS " << m_mesh.get_nx() << " " << m_mesh.get_ny() << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << m_mesh.get_dx() << " " << m_mesh.get_dy() << " 1\n"; 
    out << "POINT_DATA " << m_mesh.get_nx() * m_mesh.get_ny() << "\n";

    // Écriture des champs
    for (size_t i = 0; i < field_names.size(); ++i) {
        out << "SCALARS " << field_names[i] << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        for (const auto& val : *data_fields[i]) {
            out << val << "\n";
        }
    }
    out.close();
}

void ResultsExporter::exportStructuredVectorVTK(const std::string& filename, 
                                                const std::vector<std::string>& field_names, 
                                                const std::vector<double>& combined_data) {

    // Vérification et création du dossier
    ensure_directory_exists(filename);

    // 1. Vérification de sécurité : le vecteur doit contenir des paires (Px, Py)
    if (combined_data.size() % 2 != 0) {
        Logger::error("Le vecteur de données combinées doit contenir un nombre pair d'éléments (Px, Py).");
        return;
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }
    


    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro Vector Field\nASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << m_mesh.get_nx() << " " << m_mesh.get_ny() << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << m_mesh.get_dx() << " " << m_mesh.get_dy() << " 1\n"; 
    
    // Le nombre de points est la taille totale divisée par 2
    size_t num_points = combined_data.size() / 2;
    out << "POINT_DATA " << num_points << "\n";

    // Écriture du bloc unique de vecteurs
    out << "VECTORS " << field_names[0] << " double\n";
    for (size_t j = 0; j < combined_data.size(); j += 2) {
        out << combined_data[j] << " " << combined_data[j+1] << " 0.0\n";
    }
    out.close();
}

void ResultsExporter::exportMultiPhysicsVTK(const std::string& filename, 
                                            const std::vector<std::string>& scalar_names, 
                                            const std::vector<const std::vector<double>*>& scalar_fields,
                                            const std::vector<std::string>& vector_names, 
                                            const std::vector<const std::vector<double>*>& vector_fields) {

    ensure_directory_exists(filename);

    std::ofstream out(filename);
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    // 1. Un seul en-tête global
    out << "# vtk DataFile Version 3.0\nPhaseFieldFerro Multi-Physics\nASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << m_mesh.get_nx() << " " << m_mesh.get_ny() << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << m_mesh.get_dx() << " " << m_mesh.get_dy() << " 1\n"; 
    
    int num_points = m_mesh.get_nx() * m_mesh.get_ny();
    out << "POINT_DATA " << num_points << "\n";

    // 2. Écriture de TOUS les champs scalaires (Ta logique exacte)
    for (size_t i = 0; i < scalar_names.size(); ++i) {
        out << "SCALARS " << scalar_names[i] << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        for (const auto& val : *scalar_fields[i]) {
            out << val << "\n";
        }
    }

    // 3. Écriture de TOUS les champs vectoriels (Ta logique exacte, mais en boucle)
    for (size_t i = 0; i < vector_names.size(); ++i) {
        const std::vector<double>& combined_data = *vector_fields[i];
        
        if (combined_data.size() != static_cast<size_t>(num_points * 2)) {
            Logger::error("Le vecteur " + vector_names[i] + " n'a pas la bonne taille.");
            continue;
        }

        out << "VECTORS " << vector_names[i] << " double\n";
        for (size_t j = 0; j < combined_data.size(); j += 2) {
            out << combined_data[j] << " " << combined_data[j+1] << " 0.0\n";
        }
    }

    out.close();
}

void ResultsExporter::exportConvergenceData(const std::string& filename, 
                                            const std::vector<std::string>& column_names, 
                                            const std::vector<std::vector<double>>& data_columns) {
    // Vérification et création du dossier
    ensure_directory_exists(filename);

    std::ofstream out(filename);
    if (!out.is_open()) {
        Logger::error("Impossible d'ouvrir le fichier : " + filename);
        return;
    }

    // Écriture de l'en-tête
    for (size_t i = 0; i < column_names.size(); ++i) {
        out << column_names[i];
        if (i < column_names.size() - 1) out << "\t";
    }
    out << "\n";

    // Écriture des données ligne par ligne
    size_t num_rows = data_columns[0].size();
    for (size_t row = 0; row < num_rows; ++row) {
        for (size_t col = 0; col < data_columns.size(); ++col) {
            out << data_columns[col][row];
            if (col < data_columns.size() - 1) out << "\t";
        }
        out << "\n";
    }
    out.close();
}