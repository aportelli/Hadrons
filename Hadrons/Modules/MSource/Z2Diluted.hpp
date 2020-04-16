#ifndef Hadrons_MSource_Z2Diluted_hpp_
#define Hadrons_MSource_Z2Diluted_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Z2Diluted                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class Z2DilutedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Z2DilutedPar,
                                    std::string, noise);
};

template <typename FImpl>
class TZ2Diluted: public Module<Z2DilutedPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TZ2Diluted(const std::string name);
    // destructor
    virtual ~TZ2Diluted(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Z2Diluted, TZ2Diluted<FIMPL>, MSource);

/******************************************************************************
 *                 TZ2Diluted implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TZ2Diluted<FImpl>::TZ2Diluted(const std::string name)
: Module<Z2DilutedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TZ2Diluted<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().noise};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TZ2Diluted<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZ2Diluted<FImpl>::setup(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    envCreate(std::vector<PropagatorField>, getName(), 1, 
              noise.propSize(), envGetGrid(PropagatorField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TZ2Diluted<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto &src   = envGet(std::vector<PropagatorField>, getName());

    for (unsigned int i = 0; i < noise.propSize(); ++i)
    {
        src[i] = noise.getProp(i);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Z2Diluted_hpp_
