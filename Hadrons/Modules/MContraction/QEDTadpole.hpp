#ifndef Hadrons_MContraction_QEDTadpole_hpp_
#define Hadrons_MContraction_QEDTadpole_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               QEDTadpole                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class QEDTadpolePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDTadpolePar,
                                    std::string, loop_x,
                                    std::string, loop_y,
                                    std::string, loop_z,
                                    std::string, loop_t,
                                    std::string, emField);
};

template <typename FImpl, typename VType>
class TQEDTadpole: public Module<QEDTadpolePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);

    typedef TEmFieldGenerator<VType> EmGen;
    typedef typename EmGen::GaugeField EmField;

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
 *                     TQEDTadpole implementation                             *
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
    std::vector<std::string> in = {par().loop_x, par().loop_y, par().loop_z, par().loop_t, par().emField};
    
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
    envCreateLat(EmField, getName());
    envTmp(FFT,            "fft",         1, env().getGrid());
    envTmp(EmField,        "Gk",          1, envGetGrid(EmField));
    envTmp(LatticeComplex, "tmpcomplex",  1, envGetGrid(LatticeComplex));
    envTmp(LatticeComplex, "tmpcomplex2", 1, envGetGrid(LatticeComplex));   
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

    // Fetch env variables
    const ComplexField& qx = envGet(ComplexField, par().loop_x);
    const ComplexField& qy = envGet(ComplexField, par().loop_y);
    const ComplexField& qz = envGet(ComplexField, par().loop_z);
    const ComplexField& qt = envGet(ComplexField, par().loop_t);
    const ComplexField* q[4] = {&qx, &qy, &qz, &qt};

    const EmField& Gx = envGet(EmField, par().emField);
    envGetTmp(EmField, Gk);
    envGetTmp(LatticeComplex, tmpcomplex);
    envGetTmp(LatticeComplex, tmpcomplex2);

    LOG(Message) << "Using emField '" << par().emField << "'" << std::endl;

    // Convert emField to momentum-space
    fft.FFT_all_dim(Gk, Gx, FFT::forward);

    // Feynman gauge: emField propagator is delta^{mu, nu}/k^2.
    // Therefore we only need to compute cases where mu=nu.
    EmField& out = envGet(EmField, getName());
    for (int mu=0; mu<Nd; mu++)
    {
        // Convolve with emField propagator by multiplying in momentum-space
        fft.FFT_all_dim(tmpcomplex2, *q[mu], FFT::forward);
        const auto& em_val = peekLorentz(Gk,mu);
        tmpcomplex = tmpcomplex2 * (em_val * em_val);
        fft.FFT_all_dim(tmpcomplex2,tmpcomplex,FFT::backward);
        tmpcomplex = tmpcomplex2;

        pokeLorentz(out, tmpcomplex, mu);
    }   
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDTadpole_hpp_
