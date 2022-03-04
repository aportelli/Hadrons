#ifndef Hadrons_MIO_LoadInterlacedDistillationNoise_hpp_
#define Hadrons_MIO_LoadInterlacedDistillationNoise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         LoadInterlacedDistillationNoise                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadInterlacedDistillationNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadInterlacedDistillationNoisePar,
                                    unsigned int, ti,
                                    unsigned int, li,
                                    unsigned int, si,
                                    unsigned int, nNoise,
                                    std::string, lapEigenPack,
                                    std::string, fileName,);
};

template <typename FImpl>
class TLoadInterlacedDistillationNoise: public Module<LoadInterlacedDistillationNoisePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename DistillationNoise<FImpl>::Index Index;
public:
    // constructor
    TLoadInterlacedDistillationNoise(const std::string name);
    // destructor
    virtual ~TLoadInterlacedDistillationNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadInterlacedDistillationNoise, TLoadInterlacedDistillationNoise<FIMPL>, MIO);

/******************************************************************************
 *                 TLoadInterlacedDistillationNoise implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadInterlacedDistillationNoise<FImpl>::TLoadInterlacedDistillationNoise(const std::string name)
: Module<LoadInterlacedDistillationNoisePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadInterlacedDistillationNoise<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().lapEigenPack};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadInterlacedDistillationNoise<FImpl>::getReference(void)
{
    std::vector<std::string> ref = {par().lapEigenPack};

    return ref;
}

template <typename FImpl>
std::vector<std::string> TLoadInterlacedDistillationNoise<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadInterlacedDistillationNoise<FImpl>::setup(void)
{
    auto          &epack = envGet(typename DistillationNoise<FImpl>::LapPack, 
                                    par().lapEigenPack);
    GridCartesian *g     = envGetGrid(FermionField);
    GridCartesian *g3d   = envGetSliceGrid(FermionField, g->Nd() - 1);

    envCreateDerived(DistillationNoise<FImpl>, InterlacedDistillationNoise<FImpl>,
                     getName(), 1, g, g3d, epack, par().ti, par().li, par().si, 
                     par().nNoise);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadInterlacedDistillationNoise<FImpl>::execute(void)
{
    LOG(Message) << "Loading interlaced distillation noise with (ti, li, si) = (" 
                 << par().ti << ", " << par().li << ", " << par().si << ")"
                 << std::endl;

    auto &noise = envGetDerived(DistillationNoise<FImpl>,
                                InterlacedDistillationNoise<FImpl>, getName());
    noise.load(par().fileName, "InterlacedDistillation", vm().getTrajectory());
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

#endif // Hadrons_MIO_LoadInterlacedDistillationNoise_hpp_
