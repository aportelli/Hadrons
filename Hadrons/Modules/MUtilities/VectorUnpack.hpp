#ifndef Hadrons_MUtilities_VectorUnpack_hpp_
#define Hadrons_MUtilities_VectorUnpack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                Utility module to unpack a vector of fields                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class VectorUnpackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VectorUnpackPar,
                                    std::string, input);
};

template <typename Field>
class TVectorUnpack: public Module<VectorUnpackPar>
{
public:
    // constructor
    TVectorUnpack(const std::string name);
    // destructor
    virtual ~TVectorUnpack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ComplexVectorUnpack, TVectorUnpack<FIMPL::ComplexField>, MUtilities);
MODULE_REGISTER_TMP(FermionVectorUnpack, TVectorUnpack<FIMPL::FermionField>, MUtilities);
MODULE_REGISTER_TMP(PropagatorVectorUnpack, TVectorUnpack<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                       TVectorUnpack implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TVectorUnpack<Field>::TVectorUnpack(const std::string name)
: Module<VectorUnpackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TVectorUnpack<Field>::getInput(void)
{
    std::vector<std::string> in = {par().input};
    
    return in;
}

template <typename Field>
std::vector<std::string> TVectorUnpack<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorUnpack<Field>::setup(void)
{
    auto         &vec  = envGet(std::vector<Field>, par().input);
    unsigned int Ls    = env().getObjectLs(par().input);
    auto         *grid = vec[0].Grid();

    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        envCreate(Field, getName() + "_" + std::to_string(i), Ls, grid);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorUnpack<Field>::execute(void)
{
    auto &vec = envGet(std::vector<Field>, par().input);

    LOG(Message) << "Unpacking vector '" << par().input << "'" << std::endl;
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        auto &veci = envGet(Field, getName() + "_" + std::to_string(i));

        veci = vec[i];
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_VectorUnpack_hpp_
