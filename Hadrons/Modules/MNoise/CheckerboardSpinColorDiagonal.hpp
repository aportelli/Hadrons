#ifndef Hadrons_MNoise_CheckerboardSpinColorDiagonal_hpp_
#define Hadrons_MNoise_CheckerboardSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         CheckerboardSpinColorDiagonal                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class CheckerboardSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CheckerboardSpinColorDiagonalPar,
                                    unsigned int, nsrc,
                                    unsigned int, nsparse);
};

template <typename FImpl>
class TCheckerboardSpinColorDiagonal: public Module<CheckerboardSpinColorDiagonalPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCheckerboardSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TCheckerboardSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CheckerboardSpinColorDiagonal, TCheckerboardSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZCheckerboardSpinColorDiagonal, TCheckerboardSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TCheckerboardSpinColorDiagonal implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCheckerboardSpinColorDiagonal<FImpl>::TCheckerboardSpinColorDiagonal(const std::string name)
: Module<CheckerboardSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCheckerboardSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TCheckerboardSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCheckerboardSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(SpinColorDiagonalNoise<FImpl>, 
                     CheckerboardNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc, par().nsparse);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCheckerboardSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, getName());
    LOG(Message) << "Generating checkerboard spin-color diagonal noise with" 
                 << " nsrc = " << par().nsrc
                 << " and nSparse = " << par().nsparse << std::endl;
    noise.generateNoise(rng4d());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_CheckerboardSpinColorDiagonal_hpp_
