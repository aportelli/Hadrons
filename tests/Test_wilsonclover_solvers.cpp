/*
 * Test_wilsonclover_solvers.cpp, part of Hadrons
 *
 * Compare Hadrons Wilson-clover fixed- and mixed-precision solver modules.
 */
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    Application application;

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1501;
    globalPar.trajCounter.step  = 10;
    globalPar.runId             = "wilsonclover_solvers";
    globalPar.scheduler.schedulerType = "naive";
    application.setPar(globalPar);

    application.createModule<MGauge::Random>("gauge");

    MSource::Point::Par sourcePar;
    sourcePar.position = "0 0 0 8";
    application.createModule<MSource::Point>("source", sourcePar);

    MAction::WilsonClover::Par actionPar;
    actionPar.gauge    = "gauge";
    actionPar.mass     = 0.3;
    actionPar.csw_r    = 1.2;
    actionPar.csw_t    = 1.2;
    actionPar.cF       = 1.0;
    actionPar.boundary = "1 1 1 -1";
    actionPar.twist    = "0 0 0 0";
    application.createModule<MAction::WilsonClover>("WilsonClover", actionPar);

    MUtilities::GaugeSinglePrecisionCast::Par gaugeCastPar;
    gaugeCastPar.field = "gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gauge_f", gaugeCastPar);

    actionPar.gauge = "gauge_f";
    application.createModule<MAction::WilsonCloverF>("WilsonCloverF", actionPar);

    MSolver::RBPrecBiCGSTAB::Par bicgstabPar;
    bicgstabPar.action       = "WilsonClover";
    bicgstabPar.maxIteration = 2000;
    bicgstabPar.residual     = 1.0e-12;
    bicgstabPar.guesser      = "";
    application.createModule<MSolver::RBPrecBiCGSTAB>("BiCGSTAB", bicgstabPar);

    MSolver::RBPrecFGMRES::Par fgmresPar;
    fgmresPar.action         = "WilsonClover";
    fgmresPar.maxIteration   = 2000;
    fgmresPar.restartLength  = 25;
    fgmresPar.residual       = 1.0e-12;
    fgmresPar.guesser        = "";
    application.createModule<MSolver::RBPrecFGMRES>("FGMRES", fgmresPar);

    MSolver::MixedPrecisionRBPrecFGMRES::Par mixedFgmresBicgstabPar;
    mixedFgmresBicgstabPar.innerAction         = "WilsonCloverF";
    mixedFgmresBicgstabPar.outerAction         = "WilsonClover";
    mixedFgmresBicgstabPar.maxInnerIteration   = 50;
    mixedFgmresBicgstabPar.maxOuterIteration   = 2000;
    mixedFgmresBicgstabPar.restartLength       = 25;
    mixedFgmresBicgstabPar.innerSolver         = "BiCGSTAB";
    mixedFgmresBicgstabPar.innerRestartLength  = 10;
    mixedFgmresBicgstabPar.innerResidual       = 1.0e-4;
    mixedFgmresBicgstabPar.residual            = 1.0e-12;
    mixedFgmresBicgstabPar.outerGuesser        = "";
    application.createModule<MSolver::MixedPrecisionRBPrecFGMRES>("MixedFGMRESBiCGSTAB", mixedFgmresBicgstabPar);

    auto mixedFgmresFgmresPar = mixedFgmresBicgstabPar;
    mixedFgmresFgmresPar.innerSolver        = "FGMRES";
    mixedFgmresFgmresPar.innerRestartLength = 10;
    application.createModule<MSolver::MixedPrecisionRBPrecFGMRES>("MixedFGMRESFGMRES", mixedFgmresFgmresPar);

    MSolver::MixedPrecisionRBPrecBiCGSTAB::Par mixedBicgstabPar;
    mixedBicgstabPar.innerAction         = "WilsonCloverF";
    mixedBicgstabPar.outerAction         = "WilsonClover";
    mixedBicgstabPar.maxInnerIteration   = 50;
    mixedBicgstabPar.maxOuterIteration   = 2000;
    mixedBicgstabPar.residual            = 1.0e-12;
    mixedBicgstabPar.innerGuesser        = "";
    mixedBicgstabPar.outerGuesser        = "";
    application.createModule<MSolver::MixedPrecisionRBPrecBiCGSTAB>("MixedBiCGSTAB", mixedBicgstabPar);

    mixedBicgstabPar.residual            = 1.0e-06;
    application.createModule<MSolver::MixedPrecisionRBPrecBiCGSTAB>("MixedBiCGSTAB_lowprec", mixedBicgstabPar);

    MFermion::GaugeProp::Par propagatorPar;
    propagatorPar.source = "source";
    propagatorPar.solver = "BiCGSTAB";
    application.createModule<MFermion::GaugeProp>("prop_bicgstab", propagatorPar);

    propagatorPar.solver = "FGMRES";
    application.createModule<MFermion::GaugeProp>("prop_fgmres", propagatorPar);

    propagatorPar.solver = "MixedFGMRESBiCGSTAB";
    application.createModule<MFermion::GaugeProp>("prop_mixed_fgmres_bicgstab_inner", propagatorPar);

    propagatorPar.solver = "MixedFGMRESFGMRES";
    application.createModule<MFermion::GaugeProp>("prop_mixed_fgmres_fgmres_inner", propagatorPar);

    propagatorPar.solver = "MixedBiCGSTAB";
    application.createModule<MFermion::GaugeProp>("prop_mixed_bicgstab", propagatorPar);

    propagatorPar.solver = "MixedBiCGSTAB_lowprec";
    application.createModule<MFermion::GaugeProp>("prop_mixed_bicgstab_lowprec", propagatorPar);

    MUtilities::NormCheckPropagator::Par normCheckPar;
    normCheckPar.output = "";
    normCheckPar.relativeTolerance = 1e-8;
    normCheckPar.absoluteTolerance = 1e-15;
    normCheckPar.failIfExceedTolerance = false;

    // Fixed-precision FGMRES versus fixed-precision BiCGSTAB.
    normCheckPar.field     = "prop_fgmres";
    normCheckPar.reference = "prop_bicgstab";
    application.createModule<MUtilities::NormCheckPropagator>("fixed_fgmres_vs_fixed_bicgstab", normCheckPar);

    // Mixed-precision BiCGSTAB versus fixed-precision BiCGSTAB.
    normCheckPar.field     = "prop_mixed_bicgstab";
    normCheckPar.reference = "prop_bicgstab";
    application.createModule<MUtilities::NormCheckPropagator>("mixed_bicgstab_vs_fixed_bicgstab", normCheckPar);

    // Mixed-precision FGMRES with BiCGSTAB inner versus fixed FGMRES.
    normCheckPar.field     = "prop_mixed_fgmres_bicgstab_inner";
    normCheckPar.reference = "prop_fgmres";
    application.createModule<MUtilities::NormCheckPropagator>("mixed_fgmres_bicgstab_inner_vs_fixed_fgmres", normCheckPar);

    // Mixed-precision FGMRES with FGMRES inner versus fixed FGMRES.
    normCheckPar.field     = "prop_mixed_fgmres_fgmres_inner";
    normCheckPar.reference = "prop_fgmres";
    application.createModule<MUtilities::NormCheckPropagator>("mixed_fgmres_fgmres_inner_vs_fixed_fgmres", normCheckPar);

    // Direct comparison of the two mixed-precision FGMRES inner solvers.
    normCheckPar.field     = "prop_mixed_fgmres_fgmres_inner";
    normCheckPar.reference = "prop_mixed_fgmres_bicgstab_inner";
    application.createModule<MUtilities::NormCheckPropagator>("mixed_fgmres_fgmres_inner_vs_bicgstab_inner", normCheckPar);

    // Both mixed-precision FGMRES variants versus mixed-precision BiCGSTAB.
    normCheckPar.field     = "prop_mixed_fgmres_bicgstab_inner";
    normCheckPar.reference = "prop_mixed_bicgstab";
    application.createModule<MUtilities::NormCheckPropagator>("mixed_fgmres_bicgstab_inner_vs_mixed_bicgstab", normCheckPar);

    normCheckPar.field     = "prop_mixed_fgmres_fgmres_inner";
    application.createModule<MUtilities::NormCheckPropagator>("mixed_fgmres_fgmres_inner_vs_mixed_bicgstab", normCheckPar);

    normCheckPar.failIfExceedTolerance = true;
    normCheckPar.field     = "prop_mixed_bicgstab_lowprec";
    application.createModule<MUtilities::NormCheckPropagator>("mixed_bicgstab_lowprec_vs_mixed_bicgstab", normCheckPar);

    application.saveParameterFile("TestWilsonCloverSolvers.xml");
    application.run();

    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
