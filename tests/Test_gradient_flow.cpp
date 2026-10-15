#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.runId                = "test";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);
    
    // gauge field
    application.createModule<MGauge::Random>("gauge");

    // ////////////////////////////////////////////////////////////////////////
    // gradient flow on gauge field ///////////////////////////////////////////
    MGradientFlow::WilsonFlow::Par gfPar; /* alternative actions: SymanzikFlow, 
                                                                  ZeuthenFlow,
                                                                  CustomFlow -> choose your own plaquette and rectangle coefficients
                                                                                                gfPar.c_plaq  gfPar.c_rect
                                          */
    gfPar.gauge = "gauge";
    gfPar.steps = 10;                     // total number of evolution steps to perform
    gfPar.step_size = 0.01;               // size of one step in flow time/a^2
    gfPar.meas_interval = 1;              // interval of steps at which to measure observables
    gfPar.output = "WilsonFlow";         
    application.createModule<MGradientFlow::WilsonFlow>("WilsonFlow",gfPar);
    // ////////////////////////////////////////////////////////////////////////

    // execution
    application.saveParameterFile("WilsonFlow.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
