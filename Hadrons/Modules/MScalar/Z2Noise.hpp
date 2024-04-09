#ifndef Hadrons_MScalar_Z2Noise_hpp_
#define Hadrons_MScalar_Z2Noise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Z2 Noise, here Z2 really means {-1, 1}                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

template <typename SImpl>
class TZ2Noise: public Module<NoPar>
{
public:
    BASIC_TYPE_ALIASES(SImpl,);
public:
    // constructor
    TZ2Noise(const std::string name);
    // destructor
    virtual ~TZ2Noise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Z2Noise, TZ2Noise<SIMPL>, MScalar);

/******************************************************************************
 *                        TZ2Noise implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TZ2Noise<SImpl>::TZ2Noise(const std::string name)
: Module<NoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TZ2Noise<SImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TZ2Noise<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TZ2Noise<SImpl>::setup(void)
{
    envCreateLat(ScalarField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TZ2Noise<SImpl>::execute(void)
{
    auto &eta = envGet(ScalarField, getName());

    bernoulli(rng4d(), eta);
    eta = 2.*eta - 1.;
    eta = real(eta);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_Z2Noise_hpp_
