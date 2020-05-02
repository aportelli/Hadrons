#ifndef Hadrons_MNPR_G1_hpp_
#define Hadrons_MNPR_G1_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         G1                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class G1Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(G1Par,
                                    std::string, Sin,
                                    std::string, Sout,
                                    std::string, pin,
                                    std::string, pout,
                                    std::string, gauge,
                                    std::string, output);
};

template <typename FImpl>
class TG1: public Module<G1Par>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
        public:
            GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                    SpinColourMatrix, fourq_scalar,
                    SpinColourMatrix, fourq_gamma5,
                    SpinColourMatrix, twoq_scalar,
                    SpinColourMatrix, twoq_gamma5);
    };
    // constructor
    TG1(const std::string name);
    // destructor
    virtual ~TG1(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(G1, TG1<FIMPL>, MNPR);

/******************************************************************************
 *                 TG1 implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TG1<FImpl>::TG1(const std::string name)
: Module<G1Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TG1<FImpl>::getInput(void)
{
    std::vector<std::string> in = { par().Sout, par().Sin, par().gauge };

    return in;
}

template <typename FImpl>
std::vector<std::string> TG1<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TG1<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for G1" << std::endl;

    envTmpLat(LatticeSpinColourMatrix, "bilinear");

    envTmpLat(LatticeComplex, "bilinear_phase");
    envTmpLat(LatticeComplex, "coordinate");

    envTmpLat(GaugeField, "dSdU");
    envTmpLat(LatticeColourMatrix, "div_field_strength");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TG1<FImpl>::execute(void)
{
    auto &Umu = envGet(GaugeField, par().gauge);

    PropagatorField Sin = envGet(PropagatorField, par().Sin);
    PropagatorField Sout = envGet(PropagatorField, par().Sout);

    envGetTmp(LatticeSpinColourMatrix, bilinear);

    std::vector<int> latt_size(env().getGrid()->FullDimensions().toVector());
    std::vector<Real> pin = strToVec<Real>(par().pin);
    std::vector<Real> pout = strToVec<Real>(par().pout);

    envGetTmp(LatticeComplex, bilinear_phase);
    envGetTmp(LatticeComplex, coordinate);

    envGetTmp(GaugeField, dSdU);
    envGetTmp(LatticeColourMatrix, div_field_strength);

    Result result;
    Gamma g5(Gamma::Algebra::Gamma5);

    //// Compute volume
    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    //// Compute phases for phasing propagators
    // bilinear_phase = exp(-i (pin - pout) \cdot x)
    bilinear_phase = Zero();
    for (int mu = 0; mu < Nd; mu++) {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;

        bilinear_phase += coordinate * (pin[mu] - pout[mu]);
    }
    Complex imag = Complex(0.0, 1.0);
    bilinear_phase = exp(-imag * bilinear_phase);

    IwasakiGaugeAction<FImpl> action(1.0);
    action.deriv(Umu, dSdU);
    // The implementation of IwasakiGaugeAction::deriv has an overall
    // normalization factor of 1 / (2 Nc) which is not present in Greg's code.
    // In the interest of matching Greg's code we undo this normalization
    // factor here
    dSdU = 2.0 * RealD(Nc) * dSdU;

    for (int mu = 0; mu < Nd; mu++) {
        Gamma gmu = Gamma::gmu[mu];

        // The computation of div_field_strength here is taken from Greg's
        // code; the exact discretization is motivated by the Lattice equations
        // of motion (see Greg's thesis for more detail)
        LatticeColourMatrix U = peekLorentz(Umu, mu);
        // From the implementation of IwasakiGaugeAction::deriv(Umu, dSdU) we
        // have div_field_strength(x) = [U_\mu(x) L_\mu(x)]_{Traceless Antihermetian}
        div_field_strength = peekLorentz(dSdU, mu);
        // div_field_strength(x) = 1/2 (U_\mu(x) L_\mu(x) + L_\mu(x - \hat{\mu}) U_\mu(x - \hat{\mu}))_TA
        U = adj(U) * div_field_strength * U;
        div_field_strength = 0.5 * (div_field_strength + Cshift(U, mu, -1));

        // This next line is not necessary since div_field_strength should
        // already be traceless and anti-hermetian
        //div_field_strength = Ta(div_field_strength);

        // Scalar G1
        bilinear = g5 * adj(Sout) * g5 * div_field_strength * gmu * Sin;
        bilinear = bilinear_phase * bilinear;
        result.twoq_scalar += (1.0 / volume) * sum(bilinear);

        bilinear = bilinear_phase * bilinear;
        result.fourq_scalar += (1.0 / volume) * sum(bilinear);

        // Pseudoscalar G1
        bilinear = g5 * adj(Sout) * g5 * div_field_strength * gmu * g5 * Sin;
        bilinear = bilinear_phase * bilinear;
        result.twoq_gamma5 += (1.0 / volume) * sum(bilinear);

        bilinear = bilinear_phase * bilinear;
        result.fourq_gamma5 += (1.0 / volume) * sum(bilinear);
    }
    saveResult(par().output, "G1", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MG1_G1_hpp_
