/*
 * Test_gamma3pt_compare.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
BEGIN_MODULE_NAMESPACE(MContractionReference)

// Reference implementation copied from origin/develop at cc1a68a2.
class ReferenceGamma3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ReferenceGamma3ptPar,
                                    std::string,  q1,
                                    std::string,  q2,
                                    std::string,  q3,
                                    std::string,  gamma,
                                    unsigned int, tSnk,
                                    std::string,  output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TReferenceGamma3pt: public Module<ReferenceGamma3ptPar>
{
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    FERM_TYPE_ALIASES(FImpl3, 3);
public:
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TReferenceGamma3pt(const std::string name);
    // destructor
    virtual ~TReferenceGamma3pt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    virtual void parseGammaString(std::vector<Gamma::Algebra> &gammaList);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ReferenceGamma3pt, ARG(TReferenceGamma3pt<FIMPL, FIMPL, FIMPL>), MContractionReference);

/******************************************************************************
 *                       TReferenceGamma3pt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TReferenceGamma3pt<FImpl1, FImpl2, FImpl3>::TReferenceGamma3pt(const std::string name)
: Module<ReferenceGamma3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TReferenceGamma3pt<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().q3};

    return in;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TReferenceGamma3pt<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TReferenceGamma3pt<FImpl1, FImpl2, FImpl3>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TReferenceGamma3pt<FImpl1, FImpl2, FImpl3>::setup(void)
{
    envTmpLat(LatticeComplex, "c");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
void TReferenceGamma3pt<FImpl1, FImpl2, FImpl3>::parseGammaString(std::vector<Gamma::Algebra> &gammaList)
{
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gamma.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList.push_back((Gamma::Algebra)i);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(par().gamma);
    }
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TReferenceGamma3pt<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing 3pt contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2 << "' and '"
                 << par().q3 << "', with " << par().gamma << " insertions."
                 << std::endl;

    // Initialise variables. q2 and q3 are normal propagators, q1 may be
    // sink smeared.
    auto                        &q1 = envGet(SlicedPropagator1, par().q1);
    auto                        &q2 = envGet(PropagatorField2, par().q2);
    auto                        &q3 = envGet(PropagatorField2, par().q3);
    Gamma                       g5(Gamma::Algebra::Gamma5);
    std::vector<Gamma::Algebra> gammaList;
    std::vector<TComplex>       buf;
    std::vector<Result>         result;
    int                         nt = env().getDim(Tp);


    parseGammaString(gammaList);
    result.resize(gammaList.size());
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].gamma = gammaList[i];
        result[i].corr.resize(nt);
    }

    // Extract relevant timeslice of sinked propagator q1, then contract &
    // sum over all spacial positions of gamma insertion.
    SitePropagator1 q1Snk = q1[par().tSnk];
    envGetTmp(LatticeComplex, c);
    for (unsigned int i = 0; i < result.size(); ++i)
    {
        Gamma gamma(gammaList[i]);
        c = trace(g5*q1Snk*adj(q2)*(g5*gamma)*q3);
        sliceSum(c, buf, Tp);
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result[i].corr[t] = TensorRemove(buf[t]);
        }
    }
    saveResult(par().output, "gamma3pt", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

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

    // Test 1: compare the two Gamma3pt implementations ////////////////////
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

    MContraction::Gamma3pt::Par comparisonPar;
    comparisonPar.q1 = "Qw0_l_sliced";
    comparisonPar.q2 = "Qw0_l";
    comparisonPar.q3 = "Qp4_l";
    comparisonPar.gamma = {"Gamma5 GammaX Gamma5"};
    comparisonPar.tSnk = 4;
    comparisonPar.output = "3pt/comparison_current";
    application.createModule<MContraction::Gamma3pt>("comparison_current", comparisonPar);

    MContractionReference::ReferenceGamma3ptPar referenceGammaPar;
    referenceGammaPar.q1 = comparisonPar.q1;
    referenceGammaPar.q2 = comparisonPar.q2;
    referenceGammaPar.q3 = comparisonPar.q3;
    referenceGammaPar.gamma = "GammaX";
    referenceGammaPar.tSnk = comparisonPar.tSnk;
    referenceGammaPar.output = "3pt/comparison_reference";
    application.createModule<MContractionReference::ReferenceGamma3pt>(
        "comparison_reference", referenceGammaPar);

    LOG(Message) << "Compare 3pt/comparison_current.h5 with "
                 << "3pt/comparison_reference.h5." << std::endl;

    // execution
    application.saveParameterFile("gamma3pt.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
