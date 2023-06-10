#ifndef Hadrons_MNoise_InterlacedDistillation_hpp_
#define Hadrons_MNoise_InterlacedDistillation_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Interlaced distillation noise module                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class InterlacedDistillationPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(InterlacedDistillationPar,
                                    unsigned int,   ti,
                                    unsigned int,   li,
                                    unsigned int,   si,
                                    unsigned int,   nNoise,
                                    std::string,    lapEigenPack,
                                    std::string,    fileName);
};

template <typename FImpl>
class TInterlacedDistillation: public Module<InterlacedDistillationPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TInterlacedDistillation(const std::string name);
    // destructor
    virtual ~TInterlacedDistillation(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(InterlacedDistillation, TInterlacedDistillation<FIMPL>, MNoise);
MODULE_REGISTER_TMP(ZInterlacedDistillation, TInterlacedDistillation<ZFIMPL>, MNoise);

/******************************************************************************
 *                 TInterlacedDistillation implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TInterlacedDistillation<FImpl>::TInterlacedDistillation(const std::string name)
: Module<InterlacedDistillationPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TInterlacedDistillation<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().lapEigenPack};
    
    return in;
}

template <typename FImpl>
DependencyMap TInterlacedDistillation<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().lapEigenPack, getName()});

    return dep;
}

template <typename FImpl>
std::vector<std::string> TInterlacedDistillation<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TInterlacedDistillation<FImpl>::setup(void)
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
void TInterlacedDistillation<FImpl>::execute(void)
{
    LOG(Message) << "Generating interlaced distillation noise with (ti, li, si) = (" 
                 << par().ti << ", " << par().li << ", " << par().si << ")"
                 << std::endl;

    GridCartesian *g     = envGetGrid(FermionField);
    auto &noise = envGetDerived(DistillationNoise<FImpl>,
                                InterlacedDistillationNoise<FImpl>, getName());

    noise.generateNoise(rngSerial());
    noise.dumpDilutionMap();
    auto hash = noise.generateHash();
    LOG(Message) << "Noise hashes : " << std::endl;
    for(auto& h: hash)
    {
        LOG(Message) << h << std::endl;
    }

    if(!par().fileName.empty())
    {
        makeFileDir(par().fileName, g);
        noise.save(par().fileName, "InterlacedDistillation", vm().getTrajectory());
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_InterlacedDistillation_hpp_
