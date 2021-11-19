/*
 * Test_stoch_distil.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
    std::vector<std::string> noises  = {"noiseTdil2", "noiseTfull"};
    std::vector<std::string> nNoises = {3, 5}; 
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
    MNoise::InterlacedDistillation::Par noisePar;
    noisePar.lapEigenPack = "lapevec";
    noisePar.nNoise = nNoises[0]; // only the noises for relative lines
    noisePar.ti = 2; // T=8, ti=2 => 2 diluted time sources [0 2 4 6] and [1 3 5 7]
    noisePar.li = 3; // nVec=6, li=3 => 3 diluted eigenmode sources [0 3], [1 4] and [2 5]
    noisePar.si = 4; // nS=4, si=4 => full spin dilution; sources [0], [1], [2] and [3]
    noisePar.fileName = "./noises/Tdil2";
    application.createModule<MNoise::InterlacedDistillation>(noises[0],noisePar);
    noisePar.nNoise = nNoises[1]; // additional noises for fixed lines
    noisePar.ti = 8; // T=8, ti=8 => full time dilution
    noisePar.fileName = "./noises/Tfull";
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
            perambPar.perambFileName = "./Peramb_" + flavour[i] + "_nvec6";
            perambPar.fullSolveFileName = ""; // only used for perambMode::saveSolve
            perambPar.fullSolve = ""; // only used for perambMode::loadSolve
            perambPar.distilNoise = noises[j];
            perambPar.timeSources = tSrcs[j]; // time slices to invert on
            perambPar.perambMode = MDistil::pMode::perambOnly; // compute perambulator from lap evecs, discard unsmeared solves
            perambPar.nVec = ""; // empty = match nVec in distilNoise
            perambPar.multiFileFullSolve = ""; // delete?
            application.createModule<MDistil::Perambulator>("Peramb_" + flavour[i] + "_" + noises[j] + "_nvec6", perambPar);
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
    std::vector<std::string> noisePairsRel; 
    // brute force list of momenta
    for (int i = 0; i < nNoises[0]; ++i)
    for (int j = 0; j < nNoises[1]; ++j)
    {
        std::stringstream noisePair;
	mom << i << " " << j;
	noisePairsRel.push_back(noisePair.str());
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
    	    mom << i << " " << j;
	    noisePairsFix.push_back(noisePair.str());
        }
    }
    // meson fields

    // rho-rho field - NB: to get M(\bar{rho},rho), this needs to be multiplied by gamma5!
    MDistil::DistilMesonField::Par mfPar;
    mfPar.outPath = "./kpi-mesonfields/stoch/rho-rho";
    mfPar.lapEigenPack = "lapevec";
    mfPar.leftNoise = "noiseTfull";
    mfPar.rightNoise = "noiseTfull";
    mfPar.noisePairs = noisePairsFix;
    mfPar.leftTimeSources = tSrcs[1]; // time sources
    mfPar.rightTimeSources = tSrcs[1]; // time sources
    mfPar.leftPeramb = ""; // rho = no perambulator
    mfPar.rightPeramb = ""; // rho = no perambulator
    mfPar.leftVectorStem = "";
    mfPar.rightVectorStem = "";
    mfPar.blockSize = 24;
    mfPar.cacheSize = 4;
    mfPar.onlyDiagonal = "true";
    mfPar.gamma = "all";
    //mfPar.deltaT = 0;
    mfPar.momenta = momenta;
    application.createModule<MDistil::DistilMesonField>("RhoRho_nvec6", mfPar);
   
    // phi-phi fields
    // NB: to get M(\bar{phi},phi), this needs to be multiplied by gamma5!
    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // pion and kaon fixed meson fields for 2pt functions etc 
        //MDistil::DistilMesonField::Par mfPar;
        mfPar.outPath = "./kpi-mesonfields/stoch/phi_" + flavour[i] + "-phi_l";
        mfPar.lapEigenPack = "lapevec";
        mfPar.leftNoise = "noiseTfull";
        mfPar.rightNoise = "noiseTfull";
        mfPar.noisePairs = noisePairsFix;
        mfPar.leftTimeSources = tSrcs[1]; // time sources
        mfPar.rightTimeSources = tSrcs[1]; // time sources
        mfPar.leftPeramb = "Peramb_" + flavour[i] + "_nvec6"; // kaon / pion 1st quark
        mfPar.rightPeramb = "Peramb_l_nvec6"; // kaon / pion 2nd quark = light
        mfPar.leftVectorStem = "";
        mfPar.rightVectorStem = "";
        mfPar.blockSize = 24;
        mfPar.cacheSize = 4;
        mfPar.onlyDiagonal = "true";
        mfPar.gamma = "all";
        //mfPar.deltaT = 0;
        mfPar.momenta = momenta;
        application.createModule<MDistil::DistilMesonField>("Phi_" + flavour[i] + "phi_l_nvec6", mfPar);
    }

    // kaon relative meson fields for 2pt functions etc 
    //MDistil::DistilMesonField::Par mfPar;
    mfPar.outPath = "./kpi-mesonfields/stoch/phi_s-phi_l_rel";
    mfPar.lapEigenPack = "lapevec";
    mfPar.leftNoise = "noiseTfull";
    mfPar.rightNoise = "noiseTfull";
    mfPar.noisePairs = noisePairsFix;
    mfPar.leftTimeSources = tSrcs[1]; // time sources
    mfPar.rightTimeSources = tSrcs[1]; // time sources
    mfPar.leftPeramb = "Peramb_s_nvec6"; // kaon / pion 1st quark
    mfPar.rightPeramb = "Peramb_l_nvec6"; // kaon / pion 2nd quark = light
    mfPar.leftVectorStem = "";
    mfPar.rightVectorStem = "";
    mfPar.blockSize = 24;
    mfPar.cacheSize = 4;
    mfPar.onlyDiagonal = "true";
    mfPar.gamma = "all";
    //mfPar.deltaT = 0;
    mfPar.momenta = momenta;
    application.createModule<MDistil::DistilMesonField>("Phi_" + flavour[i] + "phi_l_nvec6", mfPar);

    // phi_l-rho fields, sink-to-sink lines with time interlaced dilution
    //MDistil::DistilMesonField::Par mfPar;
    mfPar.outPath = "./kpi-mesonfields/stoch/rho_rel-phi_l";
    mfPar.lapEigenPack = "lapevec";
    mfPar.leftNoise = "noiseTdil2"; //left field has time interlaced dilution
    mfPar.rightNoise = "noiseTfull";
    mfPar.noisePairs = noisePairsRel;
    mfPar.leftTimeSources = tSrcs[0]; // empty -> invert on all time dilution vectors
    mfPar.rightTimeSources = tSrcs[1]; // time sources
    mfPar.leftPeramb = ""; // rho
    mfPar.rightPeramb = "Peramb_l_nvec6"; // phi_l
    mfPar.leftVectorStem = "";
    mfPar.rightVectorStem = "";
    mfPar.blockSize = 24;
    mfPar.cacheSize = 4;
    mfPar.onlyDiagonal = "true";
    mfPar.gamma = "all";
    //mfPar.deltaT = 0;
    mfPar.momenta = momenta;
    application.createModule<MDistil::DistilMesonField>("Phi_" + flavour[i] + "phi_l_nvec6", mfPar);

    // execution
    application.saveParameterFile("kpiStoch.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
  
    return EXIT_SUCCESS;
}
