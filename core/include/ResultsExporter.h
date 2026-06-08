#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "Mesh.h"

using namespace std;

class ResultsExporter
{
    private:
        const Mesh& m_mesh;
        
        void ensure_directory_exists(const std::string& path);

    public:
        ResultsExporter(const Mesh& msh);

        void exportToMesh(const std::string& filename);     

        void exportStructuredScalarVTK(const std::string& filename, 
                                    const std::vector<std::string>& field_names, 
                                    const std::vector<const std::vector<double>*>& data_fields);

        void exportStructuredVectorVTK(const std::string& filename, 
                                    const std::vector<std::string>& field_names, 
                                    const std::vector<double>& combined_data);

        void exportMultiPhysicsVTK(const std::string& filename, 
                                    const std::vector<std::string>& scalar_names,
                                    const std::vector<const std::vector<double>*>& scalar_fields,
                                    const std::vector<std::string>& vector_names,
                                    const std::vector<const std::vector<double>*>& vector_fields);

        void exportConvergenceData(const std::string& filename, 
                                    const std::vector<std::string>& column_names, 
                                    const std::vector<std::vector<double>>& data_columns);

};