/*
 * Test_gamma3pt.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2024
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

BEGIN_HADRONS_NAMESPACE

template <typename FImpl>
class TUnitPropagator: public Module<NoPar>
{
public:
    typedef typename FImpl::PropagatorField PropagatorField;

    TUnitPropagator(const std::string name);
    virtual ~TUnitPropagator(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    virtual void setup(void);
    virtual void execute(void);
};

template <typename FImpl>
TUnitPropagator<FImpl>::TUnitPropagator(const std::string name)
: Module<NoPar>(name)
{}

template <typename FImpl>
std::vector<std::string> TUnitPropagator<FImpl>::getInput(void)
{
    return {};
}

template <typename FImpl>
std::vector<std::string> TUnitPropagator<FImpl>::getOutput(void)
{
    return {getName()};
}

template <typename FImpl>
void TUnitPropagator<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

template <typename FImpl>
void TUnitPropagator<FImpl>::execute(void)
{
    auto &prop = envGet(PropagatorField, getName());

    prop = 1.;
}

END_HADRONS_NAMESPACE

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MTest)

MODULE_REGISTER(UnitPropagator, TUnitPropagator<FIMPL>, MTest);

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

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

    Application              application;

    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start             = 1500;
    globalPar.trajCounter.end               = 1520;
    globalPar.trajCounter.step              = 20;
    globalPar.runId                         = "test";
    globalPar.database.restoreSchedule      = false;
    application.setPar(globalPar);

    // Test 1: analytic Gamma3pt contraction /////////////////////////////////
    application.createModule<MTest::UnitPropagator>("unit_prop");
    application.createModule<MSink::Point0>("unit_sink");

    MSink::Smear::Par unitSinkPar;
    unitSinkPar.q = "unit_prop";
    unitSinkPar.sink = "unit_sink";
    application.createModule<MSink::Smear>("unit_prop_sliced", unitSinkPar);

    MContraction::Gamma3pt::Par unitGammaPar;
    unitGammaPar.q1 = "unit_prop_sliced";
    unitGammaPar.q2 = "unit_prop";
    unitGammaPar.q3 = "unit_prop";
    unitGammaPar.gamma = {"Identity Identity Identity"};
    unitGammaPar.tSnk = 4;
    unitGammaPar.output = "3pt/unit_propagator";
    application.createModule<MContraction::Gamma3pt>("unit_gamma3pt", unitGammaPar);

    LOG(Message) << "Unit-propagator expectation: zero-momentum correlator is "
                 << "4 * Nc * spatial volume at every time slice; nonzero "
                 << "momenta vanish." << std::endl;

    // Test 2: analytic momentum projection /////////////////////////////////
    MSource::MomentumPhase::Par planeWavePar;
    planeWavePar.src = "unit_prop";
    planeWavePar.mom = "1 0 0 0";
    application.createModule<MSource::MomentumPhase>("unit_plane_wave", planeWavePar);

    MUtilities::MPScalar::Par unitMomProjPar;
    unitMomProjPar.maxFourier = 1;
    application.createModule<MUtilities::MPScalar>("unit_phases", unitMomProjPar);

    MContraction::Gamma3pt::Par planeWaveGammaPar;
    planeWaveGammaPar.q1 = "unit_prop_sliced";
    planeWaveGammaPar.q2 = "unit_prop";
    planeWaveGammaPar.q3 = "unit_plane_wave";
    planeWaveGammaPar.gamma = {"Identity Identity Identity"};
    planeWaveGammaPar.tSnk = 4;
    planeWaveGammaPar.momProjector = "unit_phases";
    planeWaveGammaPar.output = "3pt/unit_plane_wave";
    application.createModule<MContraction::Gamma3pt>("unit_plane_gamma3pt", planeWaveGammaPar);

    LOG(Message) << "Plane-wave expectation: the correlator is nonzero only at "
                 << "momentum (1, 0, 0), where it equals 4 * Nc * spatial "
                 << "volume at every time slice." << std::endl;

    // Test 3: full Gamma3pt workflow example ///////////////////////////////
    std::vector<std::string> flavour = {"l", "s"};
    std::vector<double>      mass    = {.01, .04};

    // phases for momentum projection
    MUtilities::MPScalar::Par momProjPar;
    momProjPar.maxFourier = 2;
    application.createModule<MUtilities::MPScalar>("phases", momProjPar);

    // gauge field
    application.createModule<MGauge::Random>("gauge");

    // wall source
    MSource::Z2::Par z2Par;
    z2Par.tA = 0;
    z2Par.tB = 0;
    application.createModule<MSource::Z2>("z20", z2Par);

    // point source
    MSource::Point::Par pointPar;
    pointPar.position = "0 0 0 4";
    application.createModule<MSource::Point>("point4", pointPar);

    // sink at the origin in space
    application.createModule<MSink::Point0>("sink");

    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 8;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.twist = twist;
        application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);

        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "DWF_" + flavour[i];
        solverPar.residual     = 1.0e-3;  // High residual for test purposes only. Use 1.0e-8 or smaller for physics workflows.
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);

        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "z20";
        application.createModule<MFermion::GaugeProp>("Qw0_" + flavour[i], quarkPar);
        quarkPar.source = "point4";
        application.createModule<MFermion::GaugeProp>("Qp4_" + flavour[i], quarkPar);

        // sinked propagator
        MSink::Smear::Par snkPar;
        snkPar.q = "Qw0_" + flavour[i];
        snkPar.sink = "sink";
        application.createModule<MSink::Smear>("Qw0_" + flavour[i] + "_sliced", snkPar);
    }

    // Loop over all flavour combinations (it may not make sense physically)
    for (unsigned int i=0; i<flavour.size(); ++i)
        for (unsigned int j=0; j<flavour.size(); ++j)
            for (unsigned int k=0; k<flavour.size(); ++k)
            {
                MContraction::Gamma3pt::Par threePtPar;
                threePtPar.q1 = "Qw0_" + flavour[i] + "_sliced";
                threePtPar.q2 = "Qw0_" + flavour[j];
                threePtPar.q3 = "Qp4_" + flavour[k];
                threePtPar.gamma = {    // gamma matrices in order (sink, vertex, source)
                    "Gamma5 GammaX Gamma5",
                    "GammaZ GammaY GammaX"
                };
                threePtPar.tSnk = 4;
                threePtPar.momProjector = "phases";
                threePtPar.output = "3pt/" + flavour[i] + "_" + flavour[j] + "_" + flavour[k];
                application.createModule<MContraction::Gamma3pt>("3pt_" + flavour[i] + "_" + flavour[j] + "_" + flavour[k], threePtPar);
            }

    // execution
    application.saveParameterFile("gamma3pt.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
