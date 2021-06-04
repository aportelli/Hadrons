#ifndef Hadrons_MIO_LoadInterlacedNoise_hpp_
#define Hadrons_MIO_LoadInterlacedNoise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadInterlacedNoise                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadInterlacedNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadInterlacedNoisePar,
                                    // unsigned int, nNoise,
                                    std::string, lapEigenPack,
                                    std::string, fileName,);
};

template <typename FImpl>
class TLoadInterlacedNoise: public Module<LoadInterlacedNoisePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename DistillationNoise<FImpl>::Index Index;
public:
    // constructor
    TLoadInterlacedNoise(const std::string name);
    // destructor
    virtual ~TLoadInterlacedNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadInterlacedNoise, TLoadInterlacedNoise<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadInterlacedNoise implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadInterlacedNoise<FImpl>::TLoadInterlacedNoise(const std::string name)
: Module<LoadInterlacedNoisePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadInterlacedNoise<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().lapEigenPack};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadInterlacedNoise<FImpl>::getReference(void)
{
    std::vector<std::string> ref = {par().lapEigenPack};

    return ref;
}

template <typename FImpl>
std::vector<std::string> TLoadInterlacedNoise<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadInterlacedNoise<FImpl>::setup(void)
{
    auto          &epack = envGet(typename DistillationNoise<FImpl>::LapPack, 
                                    par().lapEigenPack);
    GridCartesian *g     = envGetGrid(FermionField);
    GridCartesian *g3d   = envGetSliceGrid(FermionField, g->Nd() - 1);

    envCreateDerived(DistillationNoise<FImpl>, InterlacedDistillationNoise<FImpl>,
                     getName(), 1, g, g3d, epack);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadInterlacedNoise<FImpl>::execute(void)
{

    auto &noise = envGetDerived(DistillationNoise<FImpl>,
                                InterlacedDistillationNoise<FImpl>, getName());
    noise.load(par().fileName, "InterlacedDistillation", vm().getTrajectory());

    LOG(Message) << "Loaded interlaced distillation noise with (ti, li, si) = (" 
                    << noise.getInterlacing(Index::t) << ", "
                    << noise.getInterlacing(Index::l) << ", "
                    << noise.getInterlacing(Index::s) << ")" << std::endl;
    noise.dumpDilutionMap();
    auto hash = noise.generateHash();
    LOG(Message) << "Noise hit hashes : " << std::endl;
    for(auto& h: hash)
    {
        LOG(Message) << h << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadInterlacedNoise_hpp_
