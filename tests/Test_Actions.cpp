/*
 * Test_hadrons_spectrum.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
                       SqlNotNull<std::string>, source,
                       SqlNotNull<std::string>, action);
};


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
    globalPar.database.applicationDb        = "";//"ActionApp.db";
    globalPar.database.resultDb             = "ActionResults.db";
    globalPar.database.restoreSchedule      = false;
    globalPar.database.restoreModules       = false;
    globalPar.database.restoreMemoryProfile = false;
    globalPar.database.makeStatDb           = true;
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    application.createModule<MSink::SMatPoint>("sinkMat", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";
    std::string twist = "0. 0. 0. 0.";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions

        MAction::Wilson::Par WactionPar;
        WactionPar.gauge = "gauge";
        WactionPar.mass  = mass[i];
        WactionPar.boundary = boundary;
        WactionPar.twist = twist;
        application.createModule<MAction::Wilson>("Wilson_" + flavour[i], WactionPar);

        MAction::WilsonClover::Par WCactionPar;
        WCactionPar.gauge = "gauge";
        WCactionPar.mass  = mass[i];
        WCactionPar.csw_r = 1;
        WCactionPar.csw_t = 1;
        WCactionPar.clover_anisotropy.isAnisotropic = false;
        WCactionPar.clover_anisotropy.t_direction = 3;
        WCactionPar.clover_anisotropy.xi_0 = 1.0;
        WCactionPar.clover_anisotropy.nu = 1.0;
        WCactionPar.boundary = boundary;
        WCactionPar.twist = twist;
        application.createModule<MAction::WilsonClover>("WilsonClover_" + flavour[i], WCactionPar);

        MAction::DWF::Par DactionPar;
        DactionPar.gauge = "gauge";
        DactionPar.Ls    = 12;
        DactionPar.M5    = 1.8;
        DactionPar.mass  = mass[i];
        DactionPar.boundary = boundary;
        DactionPar.twist = twist;
        application.createModule<MAction::DWF>("DWF_" + flavour[i], DactionPar);
        
        MAction::ScaledDWF::Par SDactionPar;
        SDactionPar.gauge = "gauge";
        SDactionPar.Ls    = 12;
        SDactionPar.M5    = 1.8;
        SDactionPar.scale = 1;
        SDactionPar.mass  = mass[i];
        SDactionPar.boundary = boundary;
        SDactionPar.twist = twist;
        application.createModule<MAction::ScaledDWF>("SDWF_" + flavour[i], SDactionPar);

        MAction::MobiusDWF::Par MDactionPar;
        MDactionPar.gauge = "gauge";
        MDactionPar.Ls    = 12;
        MDactionPar.M5    = 1.8;
        MDactionPar.mass  = mass[i];
        MDactionPar.b = 1;
        MDactionPar.c = 1;
        MDactionPar.boundary = boundary;
        MDactionPar.twist = twist;
        application.createModule<MAction::MobiusDWF>("MDWF_" + flavour[i], MDactionPar);

        //MAction::ZMobiusDWF::Par zMDactionPar;
        //zMDactionPar.gauge = "gauge";
        //zMDactionPar.Ls    = 10;
        //zMDactionPar.M5    = 1.8;
        //zMDactionPar.mass  = mass[i];
        //zMDactionPar.b = 1;
        //zMDactionPar.c = 1;
        //zMDactionPar.omega = {};
        //zMDactionPar.boundary = boundary;
        //zMDactionPar.twist = twist;
        //application.createModule<MAction::ZMobiusDWF>("zMDWF_" + flavour[i], zMDactionPar);

        // solvers

        MSolver::RBPrecCG::Par WsolverPar;
        WsolverPar.action       = "Wilson_" + flavour[i];
        WsolverPar.residual     = 1.0e-8;
        WsolverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_Wilson_" + flavour[i],
                                                    WsolverPar);

        MSolver::RBPrecCG::Par WCsolverPar;
        WCsolverPar.action       = "WilsonClover_" + flavour[i];
        WCsolverPar.residual     = 1.0e-8;
        WCsolverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_WilsonColver_" + flavour[i],
                                                    WCsolverPar);

        MSolver::RBPrecCG::Par DsolverPar;
        DsolverPar.action       = "DWF_" + flavour[i];
        DsolverPar.residual     = 1.0e-8;
        DsolverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_DWF_" + flavour[i],
                                                    DsolverPar);

        MSolver::RBPrecCG::Par SDsolverPar;
        SDsolverPar.action       = "SDWF_" + flavour[i];
        SDsolverPar.residual     = 1.0e-8;
        SDsolverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_SDWF_" + flavour[i],
                                                    SDsolverPar);  

        MSolver::RBPrecCG::Par MDsolverPar;
        MDsolverPar.action       = "MDWF_" + flavour[i];
        MDsolverPar.residual     = 1.0e-8;
        MDsolverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_MDWF_" + flavour[i],
                                                    MDsolverPar);  

        //MSolver::RBPrecCG::Par zMDsolverPar;
        //zMDsolverPar.action       = "zMDWF_" + flavour[i];
        //zMDsolverPar.residual     = 1.0e-8;
        //zMDsolverPar.maxIteration = 10000;
        //application.createModule<MSolver::RBPrecCG>("CG_zMDWF_" + flavour[i],
        //                                            zMDsolverPar);                                         
        
        // propagators
        MFermion::GaugeProp::Par WquarkPar;
        WquarkPar.solver = "CG_Wilson_" + flavour[i];
        WquarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_Wilson_" + flavour[i], WquarkPar);

        MFermion::GaugeProp::Par WCquarkPar;
        WCquarkPar.solver = "CG_WilsonColver_" + flavour[i];
        WCquarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_WilsonClover_" + flavour[i], WCquarkPar);

        MFermion::GaugeProp::Par DquarkPar;
        DquarkPar.solver = "CG_DWF_" + flavour[i];
        DquarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_DWF_" + flavour[i], DquarkPar);

        MFermion::GaugeProp::Par SDquarkPar;
        SDquarkPar.solver = "CG_SDWF_" + flavour[i];
        SDquarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_SDWF_" + flavour[i], SDquarkPar);

        MFermion::GaugeProp::Par MDquarkPar;
        MDquarkPar.solver = "CG_MDWF_" + flavour[i];
        MDquarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_MDWF_" + flavour[i], MDquarkPar);

        //MFermion::GaugeProp::Par zMDquarkPar;
        //zMDquarkPar.solver = "CG_zMDWF_" + flavour[i];
        //zMDquarkPar.source = "pt";
        //application.createModule<MFermion::GaugeProp>("Qpt_zMDWF_" + flavour[i], zMDquarkPar);

    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::Meson::Par mesPar;
        MesonEntry               mesEntry;
        
        mesPar.output   = "mesons/pt_Wilson_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qpt_Wilson_" + flavour[i];
        mesPar.q2       = "Qpt_Wilson_" + flavour[j];
        mesPar.gammas   = "(Gamma5 Gamma5)";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "pt";
        mesEntry.action = "Wilson";
        application.createModule<MContraction::Meson>("meson_pt_Wilson_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_pt_Wilson_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

        mesPar.output   = "mesons/pt_WilsonClover_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qpt_WilsonClover_" + flavour[i];
        mesPar.q2       = "Qpt_WilsonClover_" + flavour[j];
        mesPar.gammas   = "(Gamma5 Gamma5)";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "pt";
        mesEntry.action = "WilsonClover";
        application.createModule<MContraction::Meson>("meson_pt_WilsonClover_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_pt_WilsonClover_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

        mesPar.output   = "mesons/pt_DWF_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qpt_DWF_" + flavour[i];
        mesPar.q2       = "Qpt_DWF_" + flavour[j];
        mesPar.gammas   = "(Gamma5 Gamma5)";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "pt";
        mesEntry.action = "DWF";
        application.createModule<MContraction::Meson>("meson_pt_DWF_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_pt_DWF_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);
                        
        mesPar.output   = "mesons/pt_SDWF_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qpt_SDWF_" + flavour[i];
        mesPar.q2       = "Qpt_SDWF_" + flavour[j];
        mesPar.gammas   = "(Gamma5 Gamma5)";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "pt";
        mesEntry.action = "SDWF";
        application.createModule<MContraction::Meson>("meson_pt_SDWF_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_pt_SDWF_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

        mesPar.output   = "mesons/pt_MDWF_" + flavour[i] + flavour[j];
        mesPar.q1       = "Qpt_MDWF_" + flavour[i];
        mesPar.q2       = "Qpt_MDWF_" + flavour[j];
        mesPar.gammas   = "(Gamma5 Gamma5)";
        mesPar.sink     = "sink";
        mesEntry.q1     = flavour[i];
        mesEntry.q2     = flavour[j];
        mesEntry.source = "pt";
        mesEntry.action = "MDWF";
        application.createModule<MContraction::Meson>("meson_pt_MDWF_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
        application.setResultMetadata("meson_pt_MDWF_" + flavour[i] + flavour[j],
                                      "meson", mesEntry);

    //    mesPar.output   = "mesons/pt_zMDWF_" + flavour[i] + flavour[j];
    //    mesPar.q1       = "Qpt_zMDWF_" + flavour[i];
    //    mesPar.q2       = "Qpt_zMDWF_" + flavour[j];
    //    mesPar.gammas   = "(Gamma5 Gamma5)";
    //    mesPar.sink     = "sink";
    //    mesEntry.q1     = flavour[i];
    //    mesEntry.q2     = flavour[j];
    //    mesEntry.source = "pt";
    //    mesEntry.action = "zMDWF";
    //    application.createModule<MContraction::Meson>("meson_pt_zMDWF_"
    //                                                  + flavour[i] + flavour[j],
    //                                                  mesPar);
    //    application.setResultMetadata("meson_pt_zMDWF_" + flavour[i] + flavour[j],
    //                                  "meson", mesEntry);

    }
    
    // execution
    application.saveParameterFile("Actions.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
