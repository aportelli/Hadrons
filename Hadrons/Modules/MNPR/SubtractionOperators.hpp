#ifndef Hadrons_MNPR_SubtractionOperators_hpp_
#define Hadrons_MNPR_SubtractionOperators_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SubtractionOperators                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class SubtractionOperatorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SubtractionOperatorsPar,
                                    std::string, Sin,
                                    std::string, Sout,
                                    std::string, pin,
                                    std::string, pout,
                                    std::string, gauge,
                                    std::string, output);
};

template <typename FImpl>
class TSubtractionOperators: public Module<SubtractionOperatorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
        public:
            // Contains information on both the 2 and 4-quark external state
            // diagrams with the given subtraction operator
            class OperatorResult : Serializable {
                public:
                GRID_SERIALIZABLE_CLASS_MEMBERS(OperatorResult,
                        SpinColourSpinColourMatrix, fourq,
                        SpinColourMatrix, twoq);
            };
            GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                    OperatorResult, dslash_left,
                    OperatorResult, dslash_gamma5_left,
                    OperatorResult, dslash_right,
                    OperatorResult, dslash_gamma5_right,
                    OperatorResult, scalar,
                    OperatorResult, psuedoscalar);
    };
    // constructor
    TSubtractionOperators(const std::string name);
    // destructor
    virtual ~TSubtractionOperators(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // covariant derivative
    virtual void dslash(PropagatorField &in, const PropagatorField &out,
        const GaugeField &Umu);

    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SubtractionOperators, TSubtractionOperators<FIMPL>, MNPR);

/******************************************************************************
 *                 TSubtractionOperators implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSubtractionOperators<FImpl>::TSubtractionOperators(const std::string name)
: Module<SubtractionOperatorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSubtractionOperators<FImpl>::getInput(void)
{
    std::vector<std::string> in = { par().Sout, par().Sin, par().gauge };

    return in;
}

template <typename FImpl>
std::vector<std::string> TSubtractionOperators<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// Copied and modified from Julia Kettle's FourQuark/Bilinear code:
// https://github.com/Julia-Kettle/Grid/blob/1ac5498c3d5ad1dc3ba8c4ce08ec6f1f102d10e0/Hadrons/Modules/MNPR/Bilinear.hpp#L120
static inline void tensor_prod(SpinColourSpinColourMatrix &lret,
        SpinColourMatrix &a, SpinColourMatrix &b)
{
    for(int si=0; si < Ns; ++si){
    for(int sj=0; sj < Ns; ++sj){
        for (int ci=0; ci < Nc; ++ci){
        for (int cj=0; cj < Nc; ++cj){
            const ComplexD val = TensorRemove(a()(si,sj)(ci,cj));
            lret()(si,sj)(ci,cj) = val * b();
        }}
    }}
}


// Computes gamma^mu D_mu for the given input field. Currently uses the
// symmetric derivative to match Greg's code, though this could change in the
// future.
template <typename FImpl>
void TSubtractionOperators<FImpl>::dslash(PropagatorField &out, const PropagatorField &in,
        const GaugeField &Umu)
{
    assert(&out != &in);
    out = Zero();
    envGetTmp(PropagatorField, tmp);
    typename FImpl::GaugeLinkField U(env().getGrid());
    for (int mu = 0; mu < Nd; mu++) {
        // Overall formula:
        // tmp(x) = U_\mu(x) in(x + \hat{\mu}) - U_\mu^\dag(x - \hat{\mu}) in(x - \hat{\mu})
        U = peekLorentz(Umu, mu);
        tmp = FImpl::CovShiftForward(U, mu, in);
        tmp = tmp - FImpl::CovShiftBackward(U, mu, in);

        Gamma gamma_mu = Gamma::gmu[mu];
        out += gamma_mu * tmp;
    }
    out = 0.5 * out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSubtractionOperators<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for subtraction operators" << std::endl;

    envTmpLat(LatticeSpinColourMatrix, "Dslash_Sin");
    envTmpLat(LatticeSpinColourMatrix, "Dslash_Sout");

    envTmpLat(LatticeSpinColourMatrix, "bilinear");

    envTmpLat(LatticeComplex, "bilinear_phase");
    envTmpLat(LatticeComplex, "pdotxout");
    envTmpLat(LatticeComplex, "coordinate");

    envTmpLat(PropagatorField, "tmp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSubtractionOperators<FImpl>::execute(void)
{
    auto &Umu = envGet(GaugeField, par().gauge);

    PropagatorField Sin = envGet(PropagatorField, par().Sin);
    PropagatorField Sout = envGet(PropagatorField, par().Sout);

    envGetTmp(LatticeSpinColourMatrix, Dslash_Sin);
    envGetTmp(LatticeSpinColourMatrix, Dslash_Sout);

    envGetTmp(LatticeSpinColourMatrix, bilinear);

    std::vector<int> latt_size(env().getGrid()->FullDimensions().toVector());
    std::vector<Real> pin = strToVec<Real>(par().pin);
    std::vector<Real> pout = strToVec<Real>(par().pout);

    envGetTmp(LatticeComplex, bilinear_phase);
    envGetTmp(LatticeComplex, pdotxout);
    envGetTmp(LatticeComplex, coordinate);

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
    pdotxout = Zero();
    for (int mu = 0; mu < Nd; mu++) {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;

        bilinear_phase += coordinate * (pin[mu] - pout[mu]);
        pdotxout += coordinate * pout[mu];
    }
    Complex imag = Complex(0.0, 1.0);
    bilinear_phase = exp(-imag * bilinear_phase);

    //// Compute Dslash for both propagators
    dslash(Dslash_Sin, Sin, Umu);
    dslash(Dslash_Sout, Sout, Umu);

    //// Compute spectator quark for 4-quark diagrams
    bilinear = Sin * exp(-imag * pdotxout);
    SpinColourMatrix spectator = sum(bilinear);

    //// Compute results
    auto compute_result = [&] (typename Result::OperatorResult &res) {
        bilinear = bilinear_phase * bilinear;
        res.twoq = (1.0 / volume) * sum(bilinear);
        bilinear = bilinear_phase * bilinear;
        SpinColourMatrix bilinear_avg = (1.0 / volume) * sum(bilinear);
        tensor_prod(res.fourq, bilinear_avg, spectator);
    };

    // The expression we want to compute here is
    //
    //  gamma^5 * adj(D_mu Sout) gamma^5 gamma^mu Sin
    //
    // But we first anti-commute the gamma^mu to the right, which gives us
    //
    //  -gamma^5 * adj(Dslash Sout) gamma^5 Sin
    //
    // Which is what we actually compute
    //
    // dslash_left here corresponds with 'Backward covariant dslash' in Greg's
    // code.
    bilinear = g5 * adj(-Dslash_Sout) * g5 * Sin;
    compute_result(result.dslash_left);

    // Note: since gamma5^2 = 1 we can simplify
    bilinear = g5 * adj(-Dslash_Sout) * Sin;
    compute_result(result.dslash_gamma5_left);

    bilinear = g5 * adj(Sout) * g5 * Dslash_Sin;
    compute_result(result.dslash_right);

    bilinear = g5 * adj(Sout) * Dslash_Sin;
    compute_result(result.dslash_gamma5_right);

    bilinear = g5 * adj(Sout) * g5 * Sin;
    compute_result(result.scalar);

    bilinear = g5 * adj(Sout) * Sin;
    compute_result(result.psuedoscalar);

    if (par().output != "") {
        saveResult(par().output, "subtraction_ops", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNPR_SubtractionOperators_hpp_
