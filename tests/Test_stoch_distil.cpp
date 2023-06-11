/*
 * Test_stoch_distil.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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
 * Two dilution schemes, both of which use 
 *   Nvec = 6       Laplacian eigenvectors
 *   LI = 3         interlaced eigenmode mode dilution
 *   SI = 4         full dilution in spin
 * for lines connected to the source timeslice ('fixed')
 *   TI = 8         full diluton in time
 *   tSrc = "0 4"   inversions only on timeslices 0 and 4
 * for sink-to-sink lines ('relative')
 *   TI = 2         interlaced time dilution
 *   tSrc = ""      empty -> use all sources, in this case 2
 *
 * We place all two-hadron oprators on the same timeslice (no 
 * displacement) and to avoid bias, we need independent noise sources.
 * A box-diagram consists of 3 fixed and 1 relative line, and these
 * are the minimum requirements for the independent noise sources.
 *   nNoiseFix = 5  for fixed lines
 *   nNoiseRel = 3  for relative lines
 *
 * We compute the following meson fields:
 * NB: to get from M(phi,phi) to M(bar{phi},phi), multiply by gamma5
 * 2pt functions:
 *   M(phi_{l/s},phi_l)   fixed 
 *   M(rho,rho)           fixed         
 * 3pt functions: 
 *   M(rho,phi_{l/s})     fixed
 * 4pt functions:
 *   M(phi_{l/s},phi_l)   relative RHS  
 *   M(rho,phi_l)         relative LHS       
 * 
 * We use the mode 'saveSolveOnly' in the light-quark perambulator module. In 
 * addition to computing the perambulators, it also saves the unsmeared 
 * solves on disk. These can late be contracted into meson fields by 
 * passing the filename via the input parameters
 *   leftVectorStem
 *   rightVectorStem
 * to the MesonFieldFixed module. This meson field represents a local 
 * current.
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
    std::vector<std::string> noises  = {"noiseTdil2", "noiseTfull"};
    std::vector<int>         nNoises = {3, 5}; 
    std::vector<std::string> tSrcs   = {"", "0 4"}; 
    
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
    // noise policies
    MNoise::InterlacedDistillation::Par noisePar;
    noisePar.lapEigenPack = "lapevec";
    noisePar.nNoise = nNoises[0]; // only the noises for relative lines
    noisePar.ti = 2; // T=8, ti=2 => 2 diluted time sources [0 2 4 6] and [1 3 5 7]
    noisePar.li = 3; // nVec=6, li=3 => 3 diluted eigenmode sources [0 3], [1 4] and [2 5]
    noisePar.si = 4; // nS=4, si=4 => full spin dilution; sources [0], [1], [2] and [3]
    noisePar.fileName = "./kpi-stoch/noises/Tdil2";
    application.createModule<MNoise::InterlacedDistillation>(noises[0],noisePar);
    noisePar.nNoise = nNoises[1]; // additional noises for fixed lines
    noisePar.ti = 8; // T=8, ti=8 => full time dilution
    noisePar.fileName = "./kpi-stoch/noises/Tfull";
    application.createModule<MNoise::InterlacedDistillation>(noises[1],noisePar);
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
        
        for (unsigned int j = 0; j < noises.size(); ++j)
        {
            // perabmulators
            MDistil::Perambulator::Par perambPar;
            perambPar.lapEigenPack = "lapevec";
            perambPar.solver = "cg_" + flavour[i];
            perambPar.perambOutFileName = "./kpi-stoch/Peramb_" + flavour[i] + "_" + noises[j];
            perambPar.unsmSolveOutFileName = "./kpi-stoch/unsmeared_solve_" + flavour[i] + "_" + noises[j]; // only used for perambMode::saveSolveOnly
            perambPar.unsmSolve = ""; // only used for perambMode::loadSolve
            perambPar.distilNoise = noises[j];
            perambPar.timeSources = tSrcs[j]; // time slices to invert on
            perambPar.perambMode = MDistil::pMode::saveSolveOnly; // compute perambulator from lap evecs, save unsmeared solves to disk
            perambPar.nVec = ""; // empty = match nVec in distilNoise
            application.createModule<MDistil::Perambulator>("Peramb_" + flavour[i] + "_" + noises[j], perambPar);
        }
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
    // noise pairs of meson fields with relative lines - independent noises
    std::vector<std::string> noisePairsRelLeft; 
    std::vector<std::string> noisePairsRelRight; 
    // brute force list of momenta
    for (int i = 0; i < nNoises[0]; ++i)
    for (int j = 0; j < nNoises[1]; ++j)
    {
        std::stringstream noisePairL;
	noisePairL << i << " " << j;
	noisePairsRelLeft.push_back(noisePairL.str());
        std::stringstream noisePairR;
	noisePairR << j << " " << i;
	noisePairsRelRight.push_back(noisePairR.str());
    }
    // noise pairs of meson fields with fixed lines - same set, so only off-diagonal!
    std::vector<std::string> noisePairsFix; 
    // brute force list of momenta
    for (int i = 0; i < nNoises[1]; ++i)
    for (int j = 0; j < nNoises[1]; ++j)
    {
        if(i != j)
        {
            std::stringstream noisePair;
    	    noisePair << i << " " << j;
	    noisePairsFix.push_back(noisePair.str());
        }
    }
    // meson fields

    // M(rho,rho)  fixed         
    MDistil::DistilMesonFieldFixed::Par mfPar;
    mfPar.outPath = "./kpi-stoch/M-rho-rho";
    mfPar.lapEigenPack = "lapevec";
    mfPar.leftNoise = "noiseTfull";
    mfPar.rightNoise = "noiseTfull";
    mfPar.noisePairs = noisePairsFix;
    mfPar.leftTimeSources = tSrcs[1]; 
    mfPar.rightTimeSources = tSrcs[1];
    // rho = no perambulator
    mfPar.leftPeramb = ""; 
    mfPar.rightPeramb = ""; 
    mfPar.leftVectorStem = "";
    mfPar.rightVectorStem = "";
    mfPar.blockSize = 12;
    mfPar.cacheSize = 4;
    mfPar.onlyDiagonal = "true";
    mfPar.gamma = "all";
    //mfPar.deltaT = 0;
    mfPar.momenta = momenta;
    application.createModule<MDistil::DistilMesonFieldFixed>("RhoRho", mfPar);
   
    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // M(phi_{l/s},phi_l)   fixed 
        mfPar.outPath = "./kpi-stoch/M-phi_" + flavour[i] + "-phi_l-fix";
        mfPar.lapEigenPack = "lapevec";
        mfPar.leftNoise = "noiseTfull";
        mfPar.rightNoise = "noiseTfull";
        mfPar.noisePairs = noisePairsFix;
        mfPar.leftTimeSources = tSrcs[1]; 
        mfPar.rightTimeSources = tSrcs[1];
        // n1 is the noises for fixed lines
        mfPar.leftPeramb = "Peramb_" + flavour[i] + "_" + noises[1]; 
        mfPar.rightPeramb = "Peramb_l_" + noises[1]; 
        mfPar.leftVectorStem = "";
        mfPar.rightVectorStem = "";
        mfPar.blockSize = 12;
        mfPar.cacheSize = 4;
        mfPar.onlyDiagonal = "true";
        mfPar.gamma = "all";
        //mfPar.deltaT = 0;
        mfPar.momenta = momenta;
        application.createModule<MDistil::DistilMesonFieldFixed>("Phi_" + flavour[i] + "phi_l-fix", mfPar);
      
        // M(rho,phi_{l/s})     fixed
        mfPar.outPath = "./kpi-stoch/M-rho-phi_" + flavour[i] + "-fix";
        mfPar.lapEigenPack = "lapevec";
        mfPar.leftNoise = "noiseTfull";
        mfPar.rightNoise = "noiseTfull";
        mfPar.noisePairs = noisePairsFix;
        mfPar.leftTimeSources = tSrcs[1];
        mfPar.rightTimeSources = tSrcs[1];
        // rho -> no perambulator
        mfPar.leftPeramb = "";
        // n1 -> fixed
        mfPar.rightPeramb = "Peramb_" + flavour[i] + "_" + noises[1];
        mfPar.leftVectorStem = "";
        mfPar.rightVectorStem = "";
        mfPar.blockSize = 12;
        mfPar.cacheSize = 4;
        mfPar.onlyDiagonal = "false";
        mfPar.gamma = "all";
        //mfPar.deltaT = 0;
        mfPar.momenta = momenta;
        application.createModule<MDistil::DistilMesonFieldFixed>("RhoPhi_" + flavour[i] + "-fix", mfPar);

        // M(phi_{l/s},phi_l)   relative RHS  
        MDistil::DistilMesonFieldRelative::Par mfrPar;
        mfrPar.outPath = "./kpi-stoch/M-phi_" + flavour[i] + "-phi_l-rel";
        mfrPar.lapEigenPack = "lapevec";
        mfrPar.leftNoise = "noiseTfull";
        mfrPar.rightNoise = "noiseTdil2";
        mfrPar.noisePairs = noisePairsRelRight;
        mfrPar.leftTimeSources = tSrcs[1]; 
        mfrPar.rightTimeSources = tSrcs[0];
        // n0 / n1 -> noise for relative / fixed lines
        mfrPar.leftPeramb = "Peramb_" + flavour[i] + "_" + noises[1];
        mfrPar.rightPeramb = "Peramb_l_" + noises[0]; 
        //mfrPar.leftVectorStem = "";
        //mfrPar.rightVectorStem = "";
        mfrPar.blockSize = 12;
        mfrPar.cacheSize = 4;
        mfrPar.gamma = "all";
        //mfrPar.deltaT = 0;
        mfrPar.momenta = momenta;
        application.createModule<MDistil::DistilMesonFieldRelative>("Phi_" + flavour[i] + "phi_l-relative", mfrPar);
    }

    // M(rho,phi_l)         relative LHS       
    MDistil::DistilMesonFieldRelative::Par mfrPar;
    mfrPar.outPath = "./kpi-stoch/M-rho-phi_l-rel";
    mfrPar.lapEigenPack = "lapevec";
    mfrPar.leftNoise = "noiseTdil2";
    mfrPar.rightNoise = "noiseTfull";
    mfrPar.noisePairs = noisePairsRelLeft;
    mfrPar.leftTimeSources = tSrcs[0]; 
    mfrPar.rightTimeSources = tSrcs[1]; 
    mfrPar.leftPeramb = ""; 
    mfrPar.rightPeramb = "Peramb_l_" + noises[1]; 
    //mfrPar.leftVectorStem = "";
    //mfrPar.rightVectorStem = "";
    mfrPar.blockSize = 12;
    mfrPar.cacheSize = 4;
    mfrPar.gamma = "all";
    //mfrPar.deltaT = 0;
    mfrPar.momenta = momenta;
    application.createModule<MDistil::DistilMesonFieldRelative>("Rhophi_l", mfrPar);

    // execution
    application.saveParameterFile("kpiStoch.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
  
    return EXIT_SUCCESS;
}
