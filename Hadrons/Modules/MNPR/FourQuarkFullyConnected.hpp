#ifndef Hadrons_MNPR_FourQuarkFullyConnected_hpp_
#define Hadrons_MNPR_FourQuarkFullyConnected_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MNPR/FourQuark.hpp>


BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MNPR)

class FourQuarkFullyConnectedPar : Serializable
{
public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FourQuarkFullyConnectedPar,
                                        std::string, qIn,
                                        std::string, qOut,
                                        std::string, pIn,
                                        std::string, pOut,
                                        std::string, gamma_basis,
                                        std::string, output);
};


template <typename FImpl>
class TFourQuarkFullyConnected : public Module<FourQuarkFullyConnectedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,)
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, gammaA,
                                        Gamma::Algebra, gammaB,
                                        std::string,  pIn,
                                        std::string,  pOut);
    };
    typedef Correlator<Metadata, SpinColourSpinColourMatrix> Result;

    TFourQuarkFullyConnected(const std::string name);
    virtual ~TFourQuarkFullyConnected(void) {};

    virtual std::vector<std::string> getInput();
    virtual std::vector<std::string> getOutput();
    virtual void tensorprod(SpinColourSpinColourMatrixField &lret, PropagatorField &a, PropagatorField &b);

protected:
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FourQuarkFullyConnected, ARG(TFourQuarkFullyConnected<FIMPL>), MNPR);

template <typename FImpl>
void TFourQuarkFullyConnected<FImpl>::tensorprod(SpinColourSpinColourMatrixField &lret, PropagatorField &a, PropagatorField &b)
{
    // Tensor product of 2 Lattice Spin Colour Matrices
    autoView(lret_v, lret, AcceleratorWrite);
    autoView(a_v, a, AcceleratorRead);
    autoView(b_v, b, AcceleratorRead);

    accelerator_for( site, lret_v.size(), grid->Nsimd(), {
        vTComplex left;
        for(int si=0; si < Ns; ++si){
        for(int sj=0; sj < Ns; ++sj){
            for (int ci=0; ci < Nc; ++ci){
            for (int cj=0; cj < Nc; ++cj){
                left()()() = a_v[site]()(si,sj)(ci,cj);
                lret_v[site]()(si,sj)(ci,cj)=left()*b_v[site]();
            }}
        }}
    });
}


template <typename FImpl>
TFourQuarkFullyConnected<FImpl>::TFourQuarkFullyConnected(const std::string name)
    : Module<FourQuarkFullyConnectedPar>(name)
{}

template <typename FImpl>
std::vector<std::string> TFourQuarkFullyConnected<FImpl>::getInput()
{
    std::vector<std::string> in = { par().qIn, par().qOut };

    return in;
}

template <typename FImpl>
std::vector<std::string> TFourQuarkFullyConnected<FImpl>::getOutput()
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
void TFourQuarkFullyConnected<FImpl>::setup()
{
    LOG(Message) << "Running setup for four-quark diagrams module"
        << std::endl;

    envTmpLat(PropagatorField, "bilinear");
    envTmpLat(PropagatorField, "bilinear_tmp");
    envTmpLat(SpinColourSpinColourMatrixField, "lret");

    envTmpLat(ComplexField, "bilinear_phase");
    envTmpLat(ComplexField, "coordinate");
}

template <typename FImpl>
void TFourQuarkFullyConnected<FImpl>::execute()
{
    LOG(Message) << "Computing contractions '" << getName()
        << "' using source propagators '" << par().qIn << "' and '" << par().qOut << "'"
        << std::endl;

    PropagatorField &qIn = envGet(PropagatorField, par().qIn);
    PropagatorField &qOut = envGet(PropagatorField, par().qOut);

    envGetTmp(PropagatorField, bilinear);
    envGetTmp(PropagatorField, bilinear_tmp);
    envGetTmp(SpinColourSpinColourMatrixField, lret);


    std::vector<Result>         result;
    Result                      r;
    std::vector<int> latt_size(env().getGrid()->FullDimensions().toVector());
    std::vector<Real> pIn = strToVec<Real>(par().pIn);
    std::vector<Real> pOut = strToVec<Real>(par().pOut);

    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);

    envGetTmp(ComplexField, bilinear_phase);
    envGetTmp(ComplexField, coordinate);

    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    LOG(Message) << "Calculating phases" << std::endl;

    bilinear_phase = Zero();
    for (int mu = 0; mu < Nd; mu++) {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2.0 * M_PI / latt_size[mu]) * coordinate;

        bilinear_phase += coordinate * (pIn[mu] - pOut[mu]);
    }
    Complex imag = Complex(0, 1.0);
    bilinear_phase = exp(-imag * bilinear_phase);

    LOG(Message) << "Done calculating phases" << std::endl;

    LOG(Message) << "Computing diagrams" << std::endl;

    r.info.pIn  = par().pIn; // Redundant to write these into every group
    r.info.pOut = par().pOut; // Redundant to write these into every group

    auto compute_diagrams = [&](Gamma gamma_A, Gamma gamma_B, bool print = true) {
        r.info.gammaA = gamma_A.g;
        r.info.gammaB = gamma_B.g;

        if (print) {
            LOG(Message) << "Computing diagrams with GammaA = "
                << gamma_A.g << ", " << "GammaB = " << gamma_B.g
                << std::endl;
        }

        // Fully connected diagram
        bilinear = bilinear_phase * (g5 * adj(qOut) * g5 * gamma_A * qIn);

        if (gamma_A.g == gamma_B.g) {
            tensorprod(lret, bilinear, bilinear);
        }
        else {
            bilinear_tmp = bilinear_phase * (g5 * adj(qOut) * g5 * gamma_B * qIn);
            tensorprod(lret, bilinear, bilinear_tmp);
        }
        r.corr.push_back( (1.0 / volume) * sum(lret) );
        result.push_back(r);
        //This is all still quite hacky - we probably want to think about the output format a little more!
        r.corr.erase(r.corr.begin());
    };

    std::string gamma_basis = par().gamma_basis;
    if (gamma_basis == "all") {
        for (Gamma gammaA: Gamma::gall) {
            for (Gamma gammaB: Gamma::gall) {
                compute_diagrams(gammaA, gammaB);
            }
        }
    }
    else if (gamma_basis == "diagonal") {
        for (Gamma g: Gamma::gall) {
            compute_diagrams(g, g);
        }
    }
    else if (gamma_basis == "diagonal_va" || gamma_basis == "diagonal_va_sp" || gamma_basis == "diagonal_va_sp_tt") {
        for (int mu = 0; mu < 4; mu++) {
            Gamma gmu = Gamma::gmu[mu];
            Gamma gmug5 = Gamma::mul[gmu.g][Gamma::Algebra::Gamma5];
            compute_diagrams(gmu, gmu);
            compute_diagrams(gmu, gmug5);
            compute_diagrams(gmug5, gmu);
            compute_diagrams(gmug5, gmug5);
        }
        if (gamma_basis == "diagonal_va_sp" || gamma_basis == "diagonal_va_sp_tt") {
            Gamma identity = Gamma(Gamma::Algebra::Identity);

            compute_diagrams(identity, identity);
            compute_diagrams(identity, g5);
            compute_diagrams(g5, identity);
            compute_diagrams(g5, g5);
        }
        if (gamma_basis == "diagonal_va_sp_tt") {
            const std::array<const Gamma, 6> gsigma = {{
                  Gamma(Gamma::Algebra::SigmaXT),      
                  Gamma(Gamma::Algebra::SigmaXY),      
                  Gamma(Gamma::Algebra::SigmaXZ),      
                  Gamma(Gamma::Algebra::SigmaYT),
                  Gamma(Gamma::Algebra::SigmaYZ),
                  Gamma(Gamma::Algebra::SigmaZT)}};

            for (Gamma gammaA: gsigma) {
                    compute_diagrams(gammaA, gammaA);
            }

            compute_diagrams(Gamma(Gamma::Algebra::SigmaXT), Gamma(Gamma::Algebra::SigmaYZ));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaXY), Gamma(Gamma::Algebra::SigmaZT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaXZ), Gamma(Gamma::Algebra::SigmaYT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaYT), Gamma(Gamma::Algebra::SigmaXZ));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaYZ), Gamma(Gamma::Algebra::SigmaXT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaZT), Gamma(Gamma::Algebra::SigmaXY));
        }
    }
    else {
        LOG(Error) << "Error: unkown gamma_basis: '" << gamma_basis << "'"
            << std::endl;
    }

    LOG(Message) << "Done computing fully-connected diagrams" << std::endl;
    saveResult(par().output, "fourquarkfullyconnected", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
