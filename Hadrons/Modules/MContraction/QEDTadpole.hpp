#ifndef Hadrons_MContraction_QEDTadpole_hpp_
#define Hadrons_MContraction_QEDTadpole_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         QEDTadpole                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class QEDTadpolePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDTadpolePar,
                                    std::string, q,
                                    std::string, photon)
};

template <typename FImpl, typename VType>
class TQEDTadpole: public Module<QEDTadpolePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);

    typedef TEmFieldGenerator<VType> EmGen;
    typedef EmGen::GaugeField EmField;

    // constructor
    TQEDTadpole(const std::string name);
    // destructor
    virtual ~TQEDTadpole(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(QEDTadpole, ARG(TQEDTadpole<FIMPL, vComplex>), MContraction);

/******************************************************************************
 *                 TQEDTadpole implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
TQEDTadpole<FImpl, VType>::TQEDTadpole(const std::string name)
: Module<QEDTadpolePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename VType>
std::vector<std::string> TQEDTadpole<FImpl, VType>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().photon};
    
    return in;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDTadpole<FImpl, VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDTadpole<FImpl, VType>::setup(void)
{
    envTmp(FFT,            "fft",         1, env().getGrid());
    envTmp(EmField,        "Ak",          1, envGetGrid(EmField));
    envTmp(LatticeComplex, "tmpcomplex",  1, envGetGrid(FermionField));
    envTmp(LatticeComplex, "tmpcomplex2", 1, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDTadpole<FImpl, VType>::execute(void)
{
    envGetTmp(FFT, fft);
    std::vector<Gamma::Algebra> Gmu = {
      Gamma::Algebra::GammaX,
      Gamma::Algebra::GammaY,
      Gamma::Algebra::GammaZ,
      Gamma::Algebra::GammaT
    };

    LOG(Message) << "Creating Tadpole field in the Feynman Gauge" << std::endl;

    // Get position-space photon
    std::cout << "Using photon '" << par().photon << "'" << std::endl;
    const EmField& Ax = envGet(EmField, par().photon);

    // Convert photon to momentum-space
    envGetTmp(EmField, Ak);
    fft.FFT_all_dim(Ak, Ax, FFT::forward);

    // Fetch env variables
    const PropagatorField& q = envGet(PropagatorField, par().q);
    envGetTmp(LatticeComplex, tmpcomplex);
    envGetTmp(LatticeComplex, tmpcomplex2);

    // Feynman gauge: photon propagator is delta^{mu, nu}/k^2.
    // Therefore we only need to compute cases where mu=nu.
    EmField& out = envGet(EmField, getName());
    for (int mu=0; mu<4; mu++)
    {
      // Contract quark loop
      tmpcomplex = trace(q*Gamma(Gmu[mu]));

      // Convolve with photon propagator by multiplying in momentum-space
      fft.FFT_all_dim(tmpcomplex2, tmpcomplex, FFT::forward);
      const auto& em_val = peekLorentz(Ak,mu);
      tmpcomplex = tmpcomplex2 * (em_val * em_val);
      fft.FFT_all_dim(tmpcomplex2,tmpcomplex,FFT::backward);

      // Remove volume factor from FFT
      tmpcomplex2 /= env().getVolume();

      pokeLorentz(out, tmpcomplex2, mu);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDTadpole_hpp_
