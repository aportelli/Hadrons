/*
 * Test_distil.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Felix Erben <felix.erben@ed.ac.uk>
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

#include <typeinfo>
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
    std::vector<std::string> flavour = {"l", "s"};
    std::vector<double>      mass    = {.01, .04};
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start             = 1500;
    globalPar.trajCounter.end               = 1520;
    globalPar.trajCounter.step              = 20;
    globalPar.runId                         = "kpi-scattering";
    globalPar.database.applicationDb        = "kpiApp.db";
    globalPar.database.resultDb             = "kpiResults.db";
    globalPar.database.restoreSchedule      = false;
    globalPar.database.restoreModules       = false;
    globalPar.database.restoreMemoryProfile = false;
    globalPar.database.makeStatDb           = false;
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Random>("gauge");
    // stout-smeared gauge field
    MGauge::StoutSmearing::Par stoutPar;
    stoutPar.gauge       = "gauge";
    stoutPar.steps       = 3;      // iterations
    stoutPar.rho         = 0.1;    // weight
    stoutPar.orthogDim   = "3";      //time-direction on a 4D grid - std::string as it can be empty for 4D smearing!
    application.createModule<MGauge::StoutSmearing>("stoutgauge", stoutPar);
    // 3D-Laplacian eigenmodes
    MDistil::LapEvec::Par lapevecPar; 
    lapevecPar.gauge = "stoutgauge";
    // Chebyshev preconditioning - this needs tuning for every new setup
    // choose odd polyOrder (-> even polynomial)
    // choose alpha < lowest eigenvector
    // choose beta > largest eigenvector
    lapevecPar.cheby.polyOrder = 11;
    lapevecPar.cheby.alpha = 1.0;
    lapevecPar.cheby.beta = 12.5;
    // nVec is relevant, the rest are tuning parameters
    lapevecPar.lanczos.nVec = 6;
    lapevecPar.lanczos.nK = 50;
    lapevecPar.lanczos.nP = 4;
    lapevecPar.lanczos.maxIt = 1000;
    lapevecPar.lanczos.resid = 1e-2;
    lapevecPar.lanczos.irlLog = 0;
    application.createModule<MDistil::LapEvec>("lapevec",lapevecPar);
    // noise policy - exact distillation
    // automatically chooses full dilution, and uses 1s as 'noise' 
    MNoise::ExactDistillation::Par noisePar;
    noisePar.lapEigenPack = "lapevec";
    //noisePar.nVec = 6; // must be smaller or equal to lapevecPar.nVec
    application.createModule<MNoise::ExactDistillation>("exact",noisePar);
    // loop over flavours, set up action, solver, perambulator 
    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // DWF actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 2; 
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        actionPar.boundary = "1 1 1 -1"; // anti-periodic bcs in time
        actionPar.twist = "0. 0. 0. 0."; // no twist
        application.createModule<MAction::DWF>("dwf_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "dwf_" + flavour[i];
        solverPar.residual     = 1e-2;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("cg_" + flavour[i],
                                                    solverPar);
        
        // perabmulators
        MDistil::Perambulator::Par perambPar;
        perambPar.lapEigenPack = "lapevec";
        perambPar.solver = "cg_" + flavour[i];
        perambPar.perambFileName = "./Peramb_" + flavour[i] + "_nvec6";
        perambPar.fullSolveFileName = ""; // only used for perambMode::saveSolve
        perambPar.fullSolve = ""; // only used for perambMode::loadSolve
        perambPar.distilNoise = "exact";
        perambPar.timeSources = "0 4"; // time slices to invert on
        perambPar.perambMode = MDistil::pMode::perambOnly; // compute perambulator from lap evecs, discard unsmeared solves
        perambPar.nVec = ""; // empty = match nVec in distilNoise
        perambPar.multiFileFullSolve = ""; // delete?
        application.createModule<MDistil::Perambulator>("Peramb_" + flavour[i] + "_nvec6", perambPar);
    }
  
    // execution
    application.saveParameterFile("kpiExact.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
  
    return EXIT_SUCCESS;
}
