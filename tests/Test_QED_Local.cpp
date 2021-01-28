/*
 * Test_QED.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
    std::vector<std::string> flavour = {"h"}; //{"l", "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {.2}; //{.01, .04, .2  , .25 , .3  };

    unsigned int  nt    = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
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

    //stochastic photon field
    MGauge::StochEm::Par photonPar;
    photonPar.gauge = PhotonR::Gauge::feynman;
    photonPar.zmScheme = PhotonR::ZmScheme::qedL;
    application.createModule<MGauge::StochEm>("ph_field", photonPar);



    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 8;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        actionPar.twist = "0. 0. 0. 0.";
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
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i],
							 quarkPar);


	//seq sources with Local vector and photon insertion
        MSource::SeqAslash::Par seqPar_V;
        seqPar_V.q         = "Qpt_" + flavour[i];
        seqPar_V.tA        = 0;
        seqPar_V.tB        = nt-1;
        seqPar_V.mom       = "0. 0. 0. 0.";
	seqPar_V.emField   = "ph_field";
        application.createModule<MSource::SeqAslash>("Qpt_" + flavour[i] 
						    + "_seq_V_ph", seqPar_V);
        
	// seq propagator with Local vector and photon insertion
        MFermion::GaugeProp::Par quarkPar_seq_V;
        quarkPar_seq_V.solver = "CG_" + flavour[i];
        quarkPar_seq_V.source = "Qpt_" + flavour[i] + "_seq_V_ph";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i] 
						+ "_seq_V_ph_" + flavour[i], 
							quarkPar_seq_V);



	//double seq sources with Local vector and photon insertion
	//(for self energy)
        MSource::SeqAslash::Par seqPar_VV;
        seqPar_VV.q         = "Qpt_" + flavour[i] + "_seq_V_ph_" 
				+ flavour[i];
        seqPar_VV.tA        = 0;
        seqPar_VV.tB        = nt-1;
        seqPar_VV.mom       = "0. 0. 0. 0.";
        seqPar_VV.emField   = "ph_field";
        application.createModule<MSource::SeqAslash>("Qpt_" + flavour[i] 
						+ "_seq_V_ph" + flavour[i] 
						+ "_seq_V_ph", seqPar_VV);
        
	//double seq propagator with Local vector and photon insertion
        MFermion::GaugeProp::Par quarkPar_seq_VV;
        quarkPar_seq_VV.solver = "CG_" + flavour[i];
        quarkPar_seq_VV.source = "Qpt_" + flavour[i] + "_seq_V_ph" 
						+ flavour[i] + "_seq_V_ph";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i] 
						+ "_seq_V_ph_" + flavour[i] 
						+ "_seq_V_ph_" + flavour[i], 
							quarkPar_seq_VV);


    }
    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        //2pt function contraction
        MContraction::Meson::Par mesPar;
        mesPar.output  = "QED/pt_" + flavour[i] + flavour[j];
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        mesPar.gammas  = "(Gamma5 Gamma5) (Gamma5 GammaTGamma5)";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);

        //photon exchange contraction
        MContraction::Meson::Par mesPar_seq_E;
        mesPar_seq_E.output  = "QED/exchange_pt_" + flavour[i] + "_V_ph_" 
				+ flavour[i] + "__" + flavour[j] + "_V_ph_"
				+ flavour[j];
        mesPar_seq_E.q1      = "Qpt_" + flavour[i] + "_seq_V_ph_" + flavour[i];
        mesPar_seq_E.q2      = "Qpt_" + flavour[j] + "_seq_V_ph_" + flavour[j];
        mesPar_seq_E.gammas  = "(Gamma5 Gamma5) (Gamma5 GammaTGamma5)";
        mesPar_seq_E.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_exchange_pt_" 
					+ flavour[i] + "_seq_V_ph_" + flavour[i] 
					+ flavour[j] + "_seq_V_ph_" + flavour[j],
                                                      mesPar_seq_E);


        //self energy contraction
        MContraction::Meson::Par mesPar_seq_S;
        mesPar_seq_S.output  = "QED/selfenergy_pt_" + flavour[i] + "_V_ph_" 
				+ flavour[i] + "_V_ph_" + flavour[i] + "__" 
				+  flavour[j];
        mesPar_seq_S.q1      = "Qpt_" + flavour[i] + "_seq_V_ph_" + flavour[i] 
				+ "_seq_V_ph_" + flavour[i];
        mesPar_seq_S.q2      = "Qpt_" + flavour[j];
        mesPar_seq_S.gammas  = "(Gamma5 Gamma5) (Gamma5 GammaTGamma5)";
        mesPar_seq_S.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_selfenergy_pt_" 
						    + flavour[i] + "_seq_V_ph_" 
						    + flavour[i] + "_seq_V_ph_" 
						    + flavour[i] + flavour[j],
                                                       mesPar_seq_S);

    }



    // execution
    application.saveParameterFile("QED.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
