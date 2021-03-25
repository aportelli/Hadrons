#ifndef Hadrons_MNPR_FourQuarkLoop_hpp_
#define Hadrons_MNPR_FourQuarkLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>


/******************************************************************************
 *                               FourQuarkLoop                                *
 ******************************************************************************/
BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MNPR)

/* This module computes loop diagrams for 4-quark operators for use in
 * non-perturbative renormalization. It can compute two diagrams, referred to
 * as "connected loop" and "disconnected loop". Diagramatically connected loop
 * refers to the following diagram:
 *
 *    Sout,pout            Sout,pout
 *        \                   /
 *         ^     _          /
 *          \  /   \      /
 *  \Gamma_A x      |    |
 *                  ^    ^
 *  \Gamma_B x      |    |
 *          /  \ _ /      \
 *         ^                \
 *        /                   \
 *    Sin,pin              Sin,pin
 *
 * Meanwhile the disconnected loop diagram is:
 *
 *     Sout,pout                       Sout,pout
 *         \                              /
 *          ^         _ _ _ _           /
 *           \      /         \       /
 *            \    |           |     |
 *   \Gamma_A  x   x \Gamma_B  ^     ^
 *            /    |           |     |
 *           /      \ _ _ _ _ /       \
 *          ^                           \
 *         /                              \
 *     Sin,pin                          Sin,pin
 *
 * In the diagrams above both 'x's are located at the same position for a given
 * diagram. At each 'x' the \Gamma_A and \Gamma_B are insertions of
 * gamma-matrix structures. There are many possiblities for the choice(s) of
 * \Gamma_A and \Gamma_B; the parameter gamma_basis controls which are used
 * (see comments below on gamma_basis).
 *
 * This diagram does *not* compute the spectator quarks (the lines on the right
 * above). The spectator quarks may be computed independently using the
 * ExternalLeg module.
 *
 * All diagrams are phased so that they carry 0 net momentum; consequently this
 * module is suitable for renormalizing 4-quark operators in the RI/SMOM scheme
 * where pin != pout.
 *
 * Outputs:
 *   There are two primary outputs of the module, both of which are written to
 *   files as vectors of SpinColourMatrix. The two outputs are referred to as
 *   bilinear and four-quark, which correspond to diagrams with 2 and 4
 *   external legs, respectively. Note that the four-quark external state
 *   diagrams should properly be SpinColourSpinColurMatrix's, but we only
 *   output a SpinColourMatrix as the amplitude of the full diagram can be
 *   obtained by taking a tensor product with the independently obtained
 *   spectator diagram.
 *
 *   The outputs gammaA and gammaB serve as a record of which insertions where
 *   taken for \Gamma_A and \Gamma_B in the diagrams above.
 *
 * Inputs:
 *   loop_type - should be either "connected" or "disconnected". Determines
 *               whether the module computes the connected or disconnected loop
 *               diagram listed above
 *
 *   Sin, Sout - Input and output propagators
 *
 *   pin, pout - The momenta corresponding to Sin and Sout, respectively
 *
 *   a2a_loop - A propagagtor field for the loop, presumably calculated using
 *              all-to-all vectors. One possible way to generate this input is
 *              using the A2ALoop module.
 *
 *   num_a2a_vectors - The number of all-to-all vectors used to generate
 *                     'a2a_loop'. This is necessary when using the module
 *                     A2ALoop since A2ALoop doesn't know the number of vectors
 *                     used, and consequently the loop calculated there is not
 *                     an average but rather a sum. We correct for this by
 *                     dividing the loop by num_a2a_vectors.
 *
 *   gamma_basis - determines the insertions \Gamma_A and \Gamma_B; see above
 *                 for more detail and see below for what options are available
 *
 *   output - struct which determines where the outputs are placed (See
 *            FourQuarkLoopOutputPar below)
 *
 */

class FourQuarkLoopOuputPar : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourQuarkLoopOuputPar,
                                    std::string, bilinear,
                                    std::string, fourquark,
                                    std::string, gammaA,
                                    std::string, gammaB);

};

class FourQuarkLoopPar : Serializable
{
public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FourQuarkLoopPar,
                                        std::string, loop_type,
                                        std::string, Sin,
                                        std::string, Sout,
                                        std::string, pin,
                                        std::string, pout,
                                        std::string, a2a_loop,
                                        int, num_a2a_vectors,
                                        std::string, gamma_basis,
                                        FourQuarkLoopOuputPar, output);
};


template <typename FImpl1, typename FImpl2>
class TFourQuarkLoop : public Module<FourQuarkLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1)
    FERM_TYPE_ALIASES(FImpl2, 2)

    TFourQuarkLoop(const std::string name);
    virtual ~TFourQuarkLoop(void) {};

    virtual std::vector<std::string> getInput();
    virtual std::vector<std::string> getOutput();

protected:
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FourQuarkLoop, ARG(TFourQuarkLoop<FIMPL, FIMPL>), MNPR);

template <typename FImpl1, typename FImpl2>
TFourQuarkLoop<FImpl1, FImpl2>::TFourQuarkLoop(const std::string name)
    : Module<FourQuarkLoopPar>(name)
{}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TFourQuarkLoop<FImpl1, FImpl2>::getInput()
{
    std::vector<std::string> in = { par().Sin, par().Sout, par().a2a_loop };

    return in;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TFourQuarkLoop<FImpl1, FImpl2>::getOutput()
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl1, typename FImpl2>
void TFourQuarkLoop<FImpl1, FImpl2>::setup()
{
    LOG(Message) << "Running setup for four-quark diagrams module"
        << std::endl;
    if (par().loop_type != "connected" && par().loop_type != "disconnected") {
        LOG(Error) << "Unkown loop type '" << par().loop_type << "'";
        assert(0);
    }

    envTmpLat(LatticeSpinColourMatrix, "loop");
    envTmpLat(LatticeSpinColourMatrix, "bilinear");

    envTmpLat(LatticeComplex, "bilinear_phase");
    envTmpLat(LatticeComplex, "coordinate");
    if (par().loop_type == "disconnected") {
        envTmpLat(LatticeComplex, "loop_trace");
    }
}

template <typename FImpl1, typename FImpl2>
void TFourQuarkLoop<FImpl1, FImpl2>::execute()
{
    LOG(Message) << "Computing contractions '" << getName()
        << "' using source propagators '" << par().Sin << "' and '" << par().Sout << "'"
        << std::endl;

    PropagatorField1 &Sin = envGet(PropagatorField1, par().Sin);
    PropagatorField2 &Sout = envGet(PropagatorField2, par().Sout);

    LatticeSpinColourMatrix &loop_raw = envGet(LatticeSpinColourMatrix, par().a2a_loop);
    envGetTmp(LatticeSpinColourMatrix, loop);
    // The loop is supposed to be an average, but the A2ALoop module doesn't
    // divide out by the number of all-to-all vectors; therefore we must divide
    // here.
    loop = loop_raw * (1.0 / par().num_a2a_vectors);

    envGetTmp(LatticeSpinColourMatrix, bilinear);

    std::vector<int> latt_size(env().getGrid()->FullDimensions().toVector());
    std::vector<Real> pin = strToVec<Real>(par().pin);
    std::vector<Real> pout = strToVec<Real>(par().pout);

    std::vector<SpinColourMatrix> fourq_result;
    std::vector<SpinColourMatrix> twoq_result;
    std::vector<Gamma::Algebra> gammaA, gammaB;
    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);

    envGetTmp(LatticeComplex, bilinear_phase);
    envGetTmp(LatticeComplex, coordinate);

    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    LOG(Message) << "Calculating phases" << std::endl;

    bilinear_phase = Zero();
    for (int mu = 0; mu < Nd; mu++) {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2.0 * M_PI / latt_size[mu]) * coordinate;

        bilinear_phase += coordinate * (pin[mu] - pout[mu]);
    }
    Complex imag = Complex(0, 1.0);
    bilinear_phase = exp(-imag * bilinear_phase);

    LOG(Message) << "Done calculating phases" << std::endl;

    LOG(Message) << "Computing diagrams" << std::endl;

    auto result_reserve = [&](int size) {
        gammaA.reserve(size);
        gammaB.reserve(size);
        fourq_result.reserve(size);
        twoq_result.reserve(size);
    };

    const bool loop_disconnected = par().loop_type == "disconnected";

    auto compute_diagrams = [&](Gamma gamma_A, Gamma gamma_B, bool print = true) {
        gammaA.push_back(gamma_A.g);
        gammaB.push_back(gamma_B.g);

        if (print) {
            LOG(Message) << "Computing diagrams with GammaA = "
                << gamma_A.g << ", " << "GammaB = " << gamma_B.g
                << std::endl;
        }

        if (loop_disconnected) {
            // Disconnected loop diagram
            envGetTmp(LatticeComplex, loop_trace);
            loop_trace = trace(loop * gamma_B);
            bilinear = g5 * adj(Sout) * g5 * gamma_A * Sin;
            bilinear = (bilinear_phase * loop_trace) * bilinear;
        }
        else {
            // Connected loop diagram
            bilinear = bilinear_phase * (g5 * adj(Sout) * g5 * gamma_A * loop * gamma_B * Sin);
        }
        SpinColourMatrix bilinear_avg = (1.0 / volume) * sum(bilinear);
        twoq_result.push_back(bilinear_avg);

         /*The quantity we actually need out of this is
         *
         * \sum_x bilinear_x \otimes spectator exp(-2i (pin - pout) \cdot x)
         *
         * Which is the appropriately phased version of the loop diagram so
         * that diagram carries zero net momentum. Using the linearity of the
         * tensor product we only need to compute
         *
         * \sum_x bilinear_x exp(-2i (pin - pout) \cdot x)
         *
         * and we can take the tensor product later.          *
         *
         * Note that at this point bilinear already has a factor of
         * bilinear_phase = exp(-i (pin - pout) \cdot x), so we only need one
         * more factor of bilinear_phase */
        bilinear = bilinear * bilinear_phase;
        bilinear_avg = (1.0 / volume) * sum(bilinear);
        fourq_result.push_back(bilinear_avg);
    };

    /* Choices of gamma_basis
     *
     * There are currently 4 choices for gamma_basis
     *
     * all - Computes every possible structure, 256 in total
     *
     * diagonal - Computes every structure where \Gamma_A = \Gamma_B, leading
     *            to 16 structures
     *
     * diagonal_va - Computes the diagrams needed to create vector/axial
     *               structure, i.e. structures of the form
     *               (qbar \Gamma^\mu q') (qbar'' \Gamma_\mu q''') where
     *               \Gamma^\mu = \gamma^\mu or \gamma^\mu \gamma^5.
     *
     * diagona_va_sp - Computes everything in 'diagonal_va' plus diagrams with
     *                 scalar/pseudoscalar strctures, i.e. structures of the
     *                 form (qbar \Gamma q') (qbar'' \Gamma q''') where \Gamma
     *                 = 1 or \gamma^5
     *
     */

    const int num_gamma = Gamma::gall.size();
    std::string gamma_basis = par().gamma_basis;
    if (gamma_basis == "all") {
        result_reserve(num_gamma * num_gamma);
        for (Gamma gammaA: Gamma::gall) {
            for (Gamma gammaB: Gamma::gall) {
                compute_diagrams(gammaA, gammaB);
            }
        }
    }
    else if (gamma_basis == "diagonal") {
        result_reserve(num_gamma);
        for (Gamma g: Gamma::gall) {
            compute_diagrams(g, g);
        }
    }
    else if (gamma_basis == "diagonal_va" || gamma_basis == "diagonal_va_sp") {
        result_reserve(4 * 4);
        for (int mu = 0; mu < 4; mu++) {
            Gamma gmu = Gamma::gmu[mu];
            Gamma gmug5 = Gamma::mul[gmu.g][Gamma::Algebra::Gamma5];
            compute_diagrams(gmu, gmu);
            compute_diagrams(gmu, gmug5);
            compute_diagrams(gmug5, gmu);
            compute_diagrams(gmug5, gmug5);
        }
        if (gamma_basis == "diagonal_va_sp") {
            result_reserve(4 * 4 + 4);
            Gamma identity = Gamma(Gamma::Algebra::Identity);

            compute_diagrams(identity, identity);
            compute_diagrams(identity, g5);
            compute_diagrams(g5, identity);
            compute_diagrams(g5, g5);
        }
    }
    else {
        LOG(Error) << "Error: unkown gamma_basis: '" << gamma_basis << "'"
            << std::endl;
    }

    LOG(Message) << "Done computing loop diagrams" << std::endl;
    saveResult(par().output.fourquark, par().loop_type + "_loop_fourq", fourq_result);
    saveResult(par().output.bilinear, par().loop_type + "_loop_twoq", twoq_result);
    saveResult(par().output.gammaA, "gammaA", gammaA);
    saveResult(par().output.gammaB, "gammaB", gammaB);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
