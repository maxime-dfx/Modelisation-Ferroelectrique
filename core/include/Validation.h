#include "Mesh.h"
#include "Logger.h"
#include "PhaseFieldFerro.h"
#include "ResultsExporter.h"

class Validation
{
    private:
        double alpha, beta, G; // Paramètres pour la validation thermodynamique
        const Datafile& datafile;
        const Mesh& msh;
        PhaseFieldFerro& physics;
        ResultsExporter& exporter;
    public:
        Validation(const Datafile& datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter);
        void validation_thermodynamique();
        void validation_electrostatique();
        void validation_mecanique();
        void validation_fracture();
        void run_all_validations();
    };