#ifndef Hadrons_MUtilities_VectorPackRef_hpp_
#define Hadrons_MUtilities_VectorPackRef_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Pack fields as a vector of pointers                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class VectorPackRefPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VectorPackRefPar,
                                    std::vector<std::string>, fields);
};

template <typename Field>
class TVectorPackRef: public Module<VectorPackRefPar>
{
public:
    // constructor
    TVectorPackRef(const std::string name);
    // destructor
    virtual ~TVectorPackRef(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(PropagatorVectorPackRef, TVectorPackRef<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                 TVectorPackRef implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TVectorPackRef<Field>::TVectorPackRef(const std::string name)
: Module<VectorPackRefPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TVectorPackRef<Field>::getInput(void)
{
    std::vector<std::string> in = par().fields;
    
    return in;
}

template <typename Field>
std::vector<std::string> TVectorPackRef<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename Field>
DependencyMap TVectorPackRef<Field>::getObjectDependencies(void)
{
    DependencyMap dep;

    for (auto &n: par().fields)
    {
        dep.insert({n, getName()});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorPackRef<Field>::setup(void)
{
    envCreate(std::vector<Field *>, getName(), 1, par().fields.size(), nullptr);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorPackRef<Field>::execute(void)
{
    LOG(Message) << "Packing " << par().fields.size() << " field pointers" << std::endl;
    LOG(Message) << "fields: " << par().fields << std::endl;
    auto &vec = envGet(std::vector<Field *>, getName());

    for (unsigned int i = 0; i < par().fields.size(); ++i)
    {
        vec[i] = &envGet(Field, par().fields[i]);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_VectorPackRef_hpp_
