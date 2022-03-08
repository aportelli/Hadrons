/*
 * Test_exact_distil.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

/**************************************************************************
 * This test shows a typical use of stochastic distillation, computing
 * all the meson fields needed to study pi-pi and K-pi scattering.
 * We assume this is run on an 8^4 grid
 *
 *   Nvec = 6            Laplacian eigenvectors
 *   exactDistillation   noise policy that enforces full dilution in spin, modes, time
 *   timeSources = ""    inversions on every available time slice (0,1,2,3,4,5,6,7)
 *
 * The only meson fields needed are:
 *   M(phi_{l/s},rho)   fixed 
 * We don't have to rely on gamma-5 hermiticity as we are computing all relevant 
 * inversions to get the 2pt, 3pt and 4pt functions.
**************************************************************************/
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
    globalPar.runId                         = "kpi-scattering-test-exact";
    globalPar.database.applicationDb        = "kpiApp.db";
    globalPar.database.resultDb             = "kpiResults.db";
    globalPar.database.restoreSchedule      = false;
    globalPar.database.restoreModules       = false;
    globalPar.database.restoreMemoryProfile = false;
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
        application.createModule<MSolver::RBPrecCG>("cg_" + flavour[i], solverPar);
        
        // perabmulators
        MDistil::Perambulator::Par perambPar;
        perambPar.lapEigenPack = "lapevec";
        perambPar.solver = "cg_" + flavour[i];
        perambPar.perambFileName = "./Peramb_" + flavour[i] + "_nvec6";
        perambPar.fullSolveFileName = ""; // only used for perambMode::saveSolve
        perambPar.fullSolve = ""; // only used for perambMode::loadSolve
        perambPar.distilNoise = "exact";
        perambPar.timeSources = ""; // empty -> invert on all time slices
        perambPar.perambMode = MDistil::pMode::perambOnly; // compute perambulator from lap evecs, discard unsmeared solves
        perambPar.nVec = ""; // empty = match nVec in distilNoise
        application.createModule<MDistil::Perambulator>("Peramb_" + flavour[i] + "_nvec6", perambPar);
    }


    // momenta
    int momSqMax=1; // P^2 max
    int maxComp = std::sqrt(momSqMax); //floor value of sqrt(P^2) = max single compnent
    std::vector<std::string> momenta; 
    // brute force list of momenta
    for (int i = -maxComp; i <= maxComp; ++i)
    for (int j = -maxComp; j <= maxComp; ++j)
    for (int k = -maxComp; k <= maxComp; ++k)
    {
       if(i*i+j*j+k*k <= momSqMax)
       {
	  std::stringstream mom;
	  mom << i << " " << j << " " << k;
	  momenta.push_back(mom.str());
       }
    }
    // meson fields

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // phi-rho fields
        MDistil::DistilMesonFieldFixed::Par mfPar;
        mfPar.outPath = "./kpi-mesonfields/phi_" + flavour[i] + "-rho";
        mfPar.lapEigenPack = "lapevec";
        mfPar.leftNoise = "exact";
        mfPar.rightNoise = "exact";
        // trivial noise policy (=identity) for exact distillation
        mfPar.noisePairs = {""};
        // empty -> invert on all time slices
        mfPar.leftTimeSources = ""; 
        mfPar.rightTimeSources = ""; 
        mfPar.leftPeramb = "Peramb_" + flavour[i] + "_nvec6";
        // rho = no perambulator
        mfPar.rightPeramb = ""; 
        // used if vectors should be read from disk instead of computed from perambulators
        mfPar.leftVectorStem = ""; 
        mfPar.rightVectorStem = "";
        // tuning parameters 
        mfPar.blockSize = 24;
        mfPar.cacheSize = 4;
        // "true" would only compute the components needed for the 2pt functions, where dilution indices are equal dt_1 = dt_2
        mfPar.onlyDiagonal = "false";
        mfPar.gamma = "all";
        // only relevant in the onlyDiagonal = "true" option, to support displaced two-hadron operators, 3pt functions etc.
        mfPar.deltaT = "0"; 
        mfPar.momenta = momenta;
        application.createModule<MDistil::DistilMesonFieldFixed>("Phi_" + flavour[i] + "Rho_nvec6", mfPar);
    }

    // execution
    application.saveParameterFile("kpiExact.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
  
    return EXIT_SUCCESS;
}
