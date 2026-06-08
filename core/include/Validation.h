#include "Mesh.h"
#include "Logger.h"
#include "PhaseFieldFerro.h"
#include "ResultsExporter.h"

class Validation
{
    private:
        double alpha, beta, G; // Paramètres pour la validation thermodynamique

    public:
        void validation_thermodynamique(const Datafile& datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter);
        void validation_electrostatique(const Datafile& datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter);
        void validation_mecanique(const Datafile& datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter);
        void validation_fracture(const Datafile& datafile, const Mesh& msh, const PhaseFieldFerro& physics, ResultsExporter& exporter);
        void run_all_validations(Datafile datafile, const Mesh& msh, PhaseFieldFerro& physics, ResultsExporter& exporter);
    };