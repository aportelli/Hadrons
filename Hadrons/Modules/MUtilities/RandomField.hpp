#ifndef Hadrons_MUtilities_RandomField_hpp_
#define Hadrons_MUtilities_RandomField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         RandomField                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class RandomFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomFieldPar,
                                    unsigned int, Ls);
};

template <typename Field>
class TRandomField: public Module<RandomFieldPar>
{
public:
    // constructor
    TRandomField(const std::string name);
    // destructor
    virtual ~TRandomField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RandomPropagator, TRandomField<FIMPL::PropagatorField>, MUtilities);
MODULE_REGISTER_TMP(RandomFermion, TRandomField<FIMPL::FermionField>, MUtilities);
MODULE_REGISTER_TMP(RandomComplex, TRandomField<FIMPL::ComplexField>, MUtilities);

/******************************************************************************
 *                 TRandomField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TRandomField<Field>::TRandomField(const std::string name)
: Module<RandomFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TRandomField<Field>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field>
std::vector<std::string> TRandomField<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TRandomField<Field>::setup(void)
{
    if (par().Ls > 1)
    {
        envCreateLat(Field, getName(), par().Ls);
    }
    else
    {
        envCreateLat(Field, getName());
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TRandomField<Field>::execute(void)
{
    LOG(Message) << "Generating random field" << std::endl;
    LOG(Message) << "Field type: " << typeName<Field>() << std::endl;
    
    auto &field = envGet(Field, getName());

    random(rng4d(), field);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_RandomField_hpp_
