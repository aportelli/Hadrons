/*
 * Test_hadrons_source.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: James Richings <james.richings@ed.ac.uk>
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

struct MesonEntry: public SqlEntry
{
    HADRONS_SQL_FIELDS(SqlNotNull<std::string>, q1, 
                       SqlNotNull<std::string>, q2,
                       SqlNotNull<std::string>, source);
};

// Info:
// This test forms a meson two point function for evey source in hadrons (ecept jacobi smear and convolution).
// The input parameters for the sources need checking.

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
    std::vector<std::string> flavour = {"l"}; //, "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {.01}; //, .04, .2  , .25 , .3  };
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start             = 1500;
    globalPar.trajCounter.end               = 1520;
    globalPar.trajCounter.step              = 20;
    globalPar.runId                         = "test";
    globalPar.database.applicationDb        = "sourceApp.db";
    globalPar.database.resultDb             = "sourceResults.db";
    globalPar.database.restoreSchedule      = false;
    globalPar.database.restoreModules       = false;
    globalPar.database.restoreMemoryProfile = false;
    globalPar.database.makeStatDb           = true;
    application.setPar(globalPar);
    
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    
    // sources

    // Point source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);

    // Momentum source
    MSource::Momentum::Par MomPar;
    MomPar.mom = "1 0 0 0";
    application.createModule<MSource::Momentum>("mom", MomPar);

    // Momentum phase
    MSource::MomentumPhase::Par MomPhPar;
    MomPhPar.src = "pt";
    MomPhPar.mom = "1 0 0 0";
    application.createModule<MSource::MomentumPhase>("momph", MomPhPar);

    // Gauss source
    MSource::Gauss::Par GaussPar;
    GaussPar.position = "0 0 0";
    GaussPar.mom = "0 0 0 0";
    GaussPar.tA = 1;
    GaussPar.tB = 1;
    GaussPar.width = 1;
    application.createModule<MSource::Gauss>("gauss", GaussPar);

    // Wall
    MSource::Wall::Par WallPar;
    WallPar.tW = 1;
    WallPar.mom = "1. 0. 0. 0.";
    application.createModule<MSource::Wall>("wall", WallPar);

    // Z2 source
    MSource::Z2::Par z2Par;
    z2Par.tA = 0;
    z2Par.tB = 0;
    application.createModule<MSource::Z2>("z2", z2Par);
    
    // Random source
    //Msource::random::Par rndPar;
    application.createModule<MSource::Random>("rnd");

    
    // missing source: Convolution, Jacobi, Z2 Dilute



    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 12;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.twist = twist;
        application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "DWF_" + flavour[i];
        solverPar.residual     = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);
        quarkPar.source = "mom";
        application.createModule<MFermion::GaugeProp>("Qmom_" + flavour[i], quarkPar);
        quarkPar.source = "momph";
        application.createModule<MFermion::GaugeProp>("Qmomph_" + flavour[i], quarkPar);
        quarkPar.source = "gauss";
        application.createModule<MFermion::GaugeProp>("Qgauss_" + flavour[i], quarkPar);
        quarkPar.source = "wall";
        application.createModule<MFermion::GaugeProp>("Qwall_" + flavour[i], quarkPar);
        quarkPar.source = "z2";
        application.createModule<MFermion::GaugeProp>("QZ2_" + flavour[i], quarkPar);
        quarkPar.source = "rnd";
        application.createModule<MFermion::GaugeProp>("Qrnd_" + flavour[i], quarkPar);


    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::Meson::Par mesPar;
        MesonEntry               mesEntry;
        
        mesPar.output   = "mesons/pt_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qpt_" + flavour[i];
        mesPar.q2       = "Qpt_" + flavour[j];
        mesPar.gammas   = "gamma5";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "pt";
        application.createModule<MContraction::Meson>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_pt_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

        mesPar.output   = "mesons/mom_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qmom_" + flavour[i];
        mesPar.q2       = "Qmom_" + flavour[j];
        mesPar.gammas   = "gamma5";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "mom";
        application.createModule<MContraction::Meson>("meson_mom_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_mom_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);


        mesPar.output   = "mesons/momph_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qmomph_" + flavour[i];
        mesPar.q2       = "Qmomph_" + flavour[j];
        mesPar.gammas   = "gamma5";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "momph";
        application.createModule<MContraction::Meson>("meson_momph_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_momph_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

        mesPar.output   = "mesons/gauss_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qgauss_" + flavour[i];
        mesPar.q2       = "Qgauss_" + flavour[j];
        mesPar.gammas   = "gamma5";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "gauss";
        application.createModule<MContraction::Meson>("meson_gauss_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_gauss_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);                              


        mesPar.output   = "mesons/wall_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qwall_" + flavour[i];
        mesPar.q2       = "Qwall_" + flavour[j];
        mesPar.gammas   = "gamma5";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "wall";
        application.createModule<MContraction::Meson>("meson_wall_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_wall_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

        mesPar.output   = "mesons/rnd_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qrnd_" + flavour[i];
        mesPar.q2       = "Qrnd_" + flavour[j];
        mesPar.gammas   = "gamma5";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "rnd";
        application.createModule<MContraction::Meson>("meson_rnd_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_rnd_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

        mesPar.output   = "mesons/Z2_" + flavour[i] + flavour[j];
        mesPar.q1       = "QZ2_" + flavour[i];
        mesPar.q2       = "QZ2_" + flavour[j];
        mesPar.gammas   = "gamma5";
        mesPar.sink     = "sink";
        mesEntry.source = "Z2";
        application.createModule<MContraction::Meson>("meson_Z2_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_Z2_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);
    }
    
    
    // execution
    application.saveParameterFile("sources.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
