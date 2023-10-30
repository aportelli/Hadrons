/*
 * Test_hadrons_meson_3pt.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
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
    globalPar.runId                = "test_conserved_bilinear";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // gauge field
    application.createModule<MGauge::Random>("gauge");

    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 1";
    std::string twist = "0. 0. 0. 0.";

    // actions
    MAction::WilsonExpClover::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.mass  = 0.01;
    actionPar.csw_r  = 1.12;
    actionPar.csw_t  = 1.12;
    actionPar.cF  = 1.0;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::WilsonExpClover>("eWC", actionPar);

    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "eWC";
    solverPar.residual     = 1.0e-15;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG", solverPar);

    // momentum source1
    MSource::Momentum::Par momentumPar1;
    momentumPar1.mom = "2 2 0 0";
    application.createModule<MSource::Momentum>("mom1", momentumPar1);

    // momentum source2
    MSource::Momentum::Par momentumPar2;
    momentumPar2.mom = "2 0 2 0";
    application.createModule<MSource::Momentum>("mom2", momentumPar2);

    // propagator1
    MFermion::GaugeProp::Par quarkPar1;
    quarkPar1.solver = "CG";
    quarkPar1.source = "mom1";
    application.createModule<MFermion::GaugeProp>("prop1", quarkPar1);

    // propagator2
    MFermion::GaugeProp::Par quarkPar2;
    quarkPar2.solver = "CG";
    quarkPar2.source = "mom2";
    application.createModule<MFermion::GaugeProp>("prop2", quarkPar2);

    // ExternalLeg1
    MNPR::ExternalLeg::Par externalLegPar1;
    externalLegPar1.qIn = "prop1";
    externalLegPar1.pIn = "2 2 0 0";
    application.createModule<MNPR::ExternalLeg>("leg1", externalLegPar1);

    // ExternalLeg2
    MNPR::ExternalLeg::Par externalLegPar2;
    externalLegPar2.qIn = "prop2";
    externalLegPar2.pIn = "2 0 2 0";
    application.createModule<MNPR::ExternalLeg>("leg2", externalLegPar2);

    // ConservedBilinear
    MNPR::ConservedBilinear::Par ConservedBilinearPar;
    ConservedBilinearPar.action = "eWC";
    ConservedBilinearPar.qIn = "prop1";
    ConservedBilinearPar.pIn = "2 2 0 0";
    ConservedBilinearPar.qOut = "prop2";
    ConservedBilinearPar.pOut = "2 0 2 0";
    application.createModule<MNPR::ConservedBilinear>("ConservedBilinear", ConservedBilinearPar);

    // execution
    application.saveParameterFile("conserved_bilinear.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
