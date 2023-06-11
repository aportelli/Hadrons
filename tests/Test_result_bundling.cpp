/*
 * Test_result_bundling.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
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
    unsigned int  nt    = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1000;
    globalPar.trajCounter.end   = 1020;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // pt source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist    = "0. 0. 0. 0.";

    // actions
    MAction::DWF::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.Ls    = 8;
    actionPar.M5    = 1.8;
    actionPar.mass  = 0.03224;
    actionPar.boundary = boundary;
    actionPar.twist = "0. 0. 0. 0.";
    application.createModule<MAction::DWF>("DWF_s", actionPar);

    
    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action       = "DWF_s";
    solverPar.residual     = 1.0e-8;
    solverPar.maxIteration = 10000;
    application.createModule<MSolver::RBPrecCG>("CG_s",
                                                solverPar);
    
    // propagators
    MFermion::GaugeProp::Par quarkPar;
    quarkPar.solver = "CG_s";
    quarkPar.source = "pt";
    application.createModule<MFermion::GaugeProp>("Qpt_s",
                            quarkPar);


    // Two Meson contractions
    MContraction::Meson::Par mesPar;
    mesPar.output  = "resultbundle_test/original_2pt_output";
    mesPar.q1      = "Qpt_s";
    mesPar.q2      = "Qpt_s";
    mesPar.gammas  = "(Gamma5 GammaX)";
    mesPar.sink    = "sink";
    application.createModule<MContraction::Meson>("meson_pt_ss_5X",
                                                    mesPar);
    MContraction::Meson::Par mesPar2;
    mesPar2.q1      = "Qpt_s";
    mesPar2.q2      = "Qpt_s";
    mesPar2.gammas  = "(Gamma5 GammaY) (Gamma5 GammaZ)";
    mesPar2.sink    = "sink";
    mesPar2.output  = "";
    application.createModule<MContraction::Meson>("meson_pt_ss_5YZ",
                                                    mesPar2);

    // A Weak Eye contraction
    MContraction::WeakNonEye3pt::Par weakNonEyePar;
    weakNonEyePar.qLeft     = "Qpt_s";
    weakNonEyePar.qBarLeft  = "Qpt_s";
    weakNonEyePar.qRight    = "Qpt_s";
    weakNonEyePar.qBarRight = "Qpt_s";
    weakNonEyePar.gammaIn  = Gamma::Algebra::GammaT;
    weakNonEyePar.gammaOut = Gamma::Algebra::Gamma5;
    weakNonEyePar.output    = "";
    application.createModule<MContraction::WeakNonEye3pt>("meson_noneye_ss",
                                                           weakNonEyePar);

    // Demonstrate the bundling of two Serializables within a group
    MIO::ResultGroup::Par resultGroupPar;
    resultGroupPar.results = std::vector<std::string>{"meson_pt_ss_5X", "meson_pt_ss_5YZ"};
    application.createModule<MIO::ResultGroup>("Group2pts",
                                               resultGroupPar);

    // Now demonstrate bundling the group with another Serializable
    MIO::WriteResultGroup::Par writeResultGroupPar;
    writeResultGroupPar.results = std::vector<std::string>{"meson_noneye_ss", "Group2pts"};
    writeResultGroupPar.output  = "resultbundle_test/output";
    application.createModule<MIO::WriteResultGroup>("WriteToFile",
                                                    writeResultGroupPar);


    // execution
    application.saveParameterFile("ResultBundling.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
