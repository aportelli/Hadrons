#ifndef Hadrons_MNoise_ExactDistillation_hpp_
#define Hadrons_MNoise_ExactDistillation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Interlaced distillation noise module                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class ExactDistillationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactDistillationPar,
                                    unsigned int, nVec,
                                    std::string, lapEigenPack);
};

template <typename FImpl>
class TExactDistillation: public Module<ExactDistillationPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TExactDistillation(const std::string name);
    // destructor
    virtual ~TExactDistillation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExactDistillation, TExactDistillation<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZExactDistillation, TExactDistillation<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TExactDistillation implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TExactDistillation<FImpl>::TExactDistillation(const std::string name)
: Module<ExactDistillationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TExactDistillation<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TExactDistillation<FImpl>::getReference(void)
{
    std::vector<std::string> ref = {par().lapEigenPack};

    return ref;
}

template <typename FImpl>
std::vector<std::string> TExactDistillation<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExactDistillation<FImpl>::setup(void)
{
    auto          &epack = envGet(typename DistillationNoise<FImpl>::LapPack, 
                                  par().lapEigenPack);
    GridCartesian *g     = envGetGrid(FermionField);
    GridCartesian *g3d   = envGetSliceGrid(FermionField, g->Nd() - 1);

    envCreateDerived(DistillationNoise<FImpl>, InterlacedDistillationNoise<FImpl>,
                     getName(), 1, g, g3d, epack, env().getDim(Tdir), par().nVec, 4, 
                     1);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExactDistillation<FImpl>::execute(void)
{
    LOG(Message) << "Generating exact distillation dummy-noise with (Nt, Nvec, Ns) = (" 
                 << env().getDim(Tdir) << ", " << par().nVec << ", 4)"
                 << std::endl;

    auto &noise = envGetDerived(DistillationNoise<FImpl>,
                                InterlacedDistillationNoise<FImpl>, getName());

    noise.exactNoisePolicy();
    noise.dumpDilutionMap();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_ExactDistillation_hpp_
