#ifndef Hadrons_MContraction_QEDSpecs_hpp_
#define Hadrons_MContraction_QEDSpecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         QEDSpecs                                 *
 ******************************************************************************/

/*
 Universal disconnected "Specs" loop subdiagram
             ___              ___
            /   \            /   \
           |     |~~~~~~~~~~|     |
            \___/   photon   \___/
            loop1            loop2

 Tr[q1 * Gamma_{mu}](x, x) * G^{mu,nu}(x, y) * Tr[q2 * Gamma_{nu}](y, y) 

*/


BEGIN_MODULE_NAMESPACE(MContraction)

class QEDSpecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDSpecsPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, emField,
                                    std::string, output)
};


template <typename FImpl, typename VType>
class TQEDSpecs: public Module<QEDSpecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);

    typedef TEmFieldGenerator<VType> EmGen;
    typedef typename EmGen::GaugeField EmField;

    // constructor
    TQEDSpecs(const std::string name);
    // destructor
    virtual ~TQEDSpecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void) override;
    virtual std::vector<std::string> getOutput(void) override;
    virtual std::vector<std::string> getOutputFiles(void) override;
    // setup
    virtual void setup(void) override;
    // execution
    virtual void execute(void) override;
};

MODULE_REGISTER_TMP(QEDSpecs, ARG(TQEDSpecs<FIMPL, vComplex>), MContraction);

/******************************************************************************
 *                 TQEDSpecs implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
TQEDSpecs<FImpl, VType>::TQEDSpecs(const std::string name)
: Module<QEDSpecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename VType>
std::vector<std::string> TQEDSpecs<FImpl, VType>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().emField};
    
    return in;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDSpecs<FImpl, VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDSpecs<FImpl, VType>::getOutputFiles(void)
{
    std::vector<std::string> output = {};
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDSpecs<FImpl, VType>::setup(void)
{
    envTmp(FFT,            "fft",         1, env().getGrid());
    envTmp(EmField,        "weight_k",    1, envGetGrid(EmField));
    envTmp(LatticeComplex, "tmpcomplex",  1, envGetGrid(FermionField));
    envTmp(LatticeComplex, "tmpcomplex2", 1, envGetGrid(FermionField));
    envTmp(LatticeReal,    "tmpreal",     1, envGetGrid(FermionField));
    envTmp(LatticeReal,    "tmpreal2",    1, envGetGrid(FermionField));
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDSpecs<FImpl, VType>::execute(void)
{
    envGetTmp(FFT, fft);
    std::vector<Gamma> Gmu = {
      Gamma(Gamma::Algebra::GammaX),
      Gamma(Gamma::Algebra::GammaY),
      Gamma(Gamma::Algebra::GammaZ),
      Gamma(Gamma::Algebra::GammaT)
    };

    LOG(Message) << "Starting Specs contraction in Feynman Gauge" << std::endl;

    // Get position-space emField
    std::cout << "Using emField '" << par().emField << "'" << std::endl;
    const EmField& weight = envGet(EmField, par().emField);

    // Convert emField to momentum-space
    envGetTmp(EmField, weight_k);
    fft.FFT_all_dim(weight_k, weight, FFT::forward);

    // Fetch env variables
    const PropagatorField& q1 = envGet(PropagatorField, par().q1);
    const PropagatorField& q2 = envGet(PropagatorField, par().q2);
    envGetTmp(LatticeComplex, tmpcomplex);
    envGetTmp(LatticeComplex, tmpcomplex2);
    envGetTmp(LatticeReal, tmpreal);
    envGetTmp(LatticeReal, tmpreal2);

    // Output variable
    RealD result = 0;

    // Feynman gauge: photon propagator is delta^{mu, nu}/k^2.
    // Therefore we only need to compute cases where mu=nu.
    for (int mu=0; mu<4; mu++)
    {
      // Contract quark loop 1
      tmpcomplex = trace(q1*Gmu[mu]);

      // Convolve with photon propagator by multiplying in momentum-space
      fft.FFT_all_dim(tmpcomplex2, tmpcomplex, FFT::forward);
      const auto& em_val = peekLorentz(weight_k,mu);
      tmpcomplex = tmpcomplex2 * (em_val * em_val);
      fft.FFT_all_dim(tmpcomplex2,tmpcomplex,FFT::backward);
      
      tmpreal = toReal(imag(tmpcomplex2));

      // Contract loop 2
      tmpcomplex = trace(q2*Gmu[mu]);
      tmpreal2   = toReal(imag(tmpcomplex));

      result += sum(tmpreal * tmpreal2);
    }

    LOG(Message) << "specs: " << result << std::endl;

    saveResult(par().output, "specs", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDSpecs_hpp_
