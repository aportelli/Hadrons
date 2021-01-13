#ifndef Hadrons_MNoise_InterlacedDilutionNoise_hpp_
#define Hadrons_MNoise_InterlacedDilutionNoise_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>
#include <Hadrons/Modules/MDistil/Distil.hpp> // TODO: for 3D grids, shodul be handled globally

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Interlaced distillation noise module                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class InterlacedDilutionNoisePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(InterlacedDilutionNoisePar,
                                    unsigned int, ti,
                                    unsigned int, li,
                                    unsigned int, si,
                                    unsigned int, nNoise,
                                    std::string, lapEigenPack);
};

template <typename FImpl>
class TInterlacedDilutionNoise: public Module<InterlacedDilutionNoisePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TInterlacedDilutionNoise(const std::string name);
    // destructor
    virtual ~TInterlacedDilutionNoise(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    // 3D grid, not great, needs to come from the environment
    std::unique_ptr<GridCartesian> g3d_;
};

MODULE_REGISTER_TMP(InterlacedDilutionNoise, TInterlacedDilutionNoise<FIMPL>, MNoise);

/******************************************************************************
 *                 TInterlacedDilutionNoise implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TInterlacedDilutionNoise<FImpl>::TInterlacedDilutionNoise(const std::string name)
: Module<InterlacedDilutionNoisePar>(name)
{
    MDistil::MakeLowerDimGrid(g3d_, envGetGrid(FermionField));
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TInterlacedDilutionNoise<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TInterlacedDilutionNoise<FImpl>::getReference(void)
{
    std::vector<std::string> ref = {par().lapEigenPack};

    return ref;
}

template <typename FImpl>
std::vector<std::string> TInterlacedDilutionNoise<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TInterlacedDilutionNoise<FImpl>::setup(void)
{
    auto &epack = envGet(typename DistillationNoise<FImpl>::LapPack, 
                         par().lapEigenPack);
    
    envCreateDerived(DistillationNoise<FImpl>, InterlacedDistillationNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), g3d_.get(), 
                     epack, par().ti, par().li, par().si, par().nNoise);

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TInterlacedDilutionNoise<FImpl>::execute(void)
{
    LOG(Message) << "Generating interlaced distillation noise with (ti, li, si) = (" 
                 << par().ti << ", " << par().li << ", " << par().si << ")"
                 << std::endl;

    auto &noise = envGetDerived(DistillationNoise<FImpl>,
                                InterlacedDistillationNoise<FImpl>, getName());

    noise.generateNoise(rngSerial());
    noise.dumpDilutionMap();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_InterlacedDilutionNoise_hpp_
