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
    double        mass    = .25;
    
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

    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    // point source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    std:: string srcName  = "pt_0";
    application.createModule<MSource::Point>(srcName, ptPar);

    // point sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // action
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = 12;
    actionPar.M5    = 1.8;
    actionPar.mass  = mass;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>("DWF", actionPar);
    
    // solver
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "DWF";
    solverPar.residual     = 1.0e-8;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG",solverPar);

    // propagator
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG";
    quarkPar.source = srcName;
    std::vector<std::string> qName = {"Qpt_0"};
    application.createModule<MFermion::GaugeProp>(qName[0], quarkPar);

    // zero flow contractions
    MContraction::Meson::Par mesPar;
    mesPar.q1     = qName[0];
    mesPar.q2     = qName[0];
    mesPar.gammas = "all";
    mesPar.sink   = "sink";
    application.createModule<MContraction::Meson>("meson_t0.00",mesPar);

    std::vector<std::string> results = {"meson_t0.00"}; // collect names of results to be written to file
    
    // ///////////////////////////////////////////////////////////////////////
    // gradient flow on propagator ///////////////////////////////////////////
    MGradientFlow::WilsonFermionFlow::Par wfPar;
    wfPar.gauge = "gauge";
    wfPar.steps = 10;                        // total number of evolution steps to perform
    wfPar.step_size = 0.01;                  // size of one step in flow time/a^2
    wfPar.meas_interval = 10;                // interval of steps at which to measure observables
    wfPar.props = qName;                     // provide a list of propagators to be flowed
    //wfPar.outProps = {"flowedquark1"};     // (optional) provide a list of names for the output propagators
    wfPar.bc = -1;                           // set boundary conditions to anti-periodic in time (bc = 1 will keep periodic)
    //wfPar.output = "GaugeFlow";            // option to output gauge flow data separately
    application.createModule<MGradientFlow::WilsonFermionFlow>("FermionFlow",wfPar);
    // ///////////////////////////////////////////////////////////////////////
    
    // positive flow contractions
    MContraction::Meson::Par mesflowPar;
    mesflowPar.q1     = qName[0]+"_t0.10";
    mesflowPar.q2     = qName[0]+"_t0.10";
    mesflowPar.gammas = "all";
    mesflowPar.sink   = "sink";
    application.createModule<MContraction::Meson>("meson_t0.10",mesflowPar);
    results.push_back("meson_t0.10");

    // save data
    MIO::WriteResultGroup::Par wPar;
    wPar.results = results;
    wPar.output = "results";
    application.createModule<MIO::WriteResultGroup>("WriteToFile",wPar);

    // execution
    application.saveParameterFile("FermionFlow.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
