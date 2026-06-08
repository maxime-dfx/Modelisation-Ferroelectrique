#pragma once
#include "PhaseFieldFerro.h"
#include "ResultsExporter.h"
#include "Datafile.h"
#include <string>

class Simulation {
private:
    PhaseFieldFerro& m_physics;
    ResultsExporter& m_exporter;
    
    int m_max_steps;
    int m_save_frequency;
    int m_current_step;
    double m_tolerance;
    std::string m_output_dir;
    std::string m_vtk_output_dir;
    std::string m_convergence_data;

    void save_current_state();

public:
    Simulation(PhaseFieldFerro& physics, ResultsExporter& exporter, const Datafile& config);
    void run();
};