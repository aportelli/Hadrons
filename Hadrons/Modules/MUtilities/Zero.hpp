#ifndef Hadrons_MUtilities_Zero_hpp_
#define Hadrons_MUtilities_Zero_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Zero                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class ZeroPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ZeroPar,
                                    unsigned int, Ls);
};


template <typename Field>
class TZero: public Module<ZeroPar>
{
public:
    // constructor
    TZero(const std::string name);
    // destructor
    virtual ~TZero(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ZeroScalarField,                     TZero<FIMPL::Field          >, MUtilities);
MODULE_REGISTER_TMP(ZeroComplexField,                    TZero<FIMPL::ComplexField   >, MUtilities);
MODULE_REGISTER_TMP(ZeroFermionField,                    TZero<FIMPL::FermionField   >, MUtilities);
MODULE_REGISTER_TMP(ZeroPropagatorField,                 TZero<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                 TZero implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TZero<Field>::TZero(const std::string name)
: Module<ZeroPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TZero<Field>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field>
std::vector<std::string> TZero<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TZero<Field>::setup(void)
{
    envCreate(Field, getName(), par().Ls, envGetGrid(Field));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TZero<Field>::execute(void)
{
    Field& out = envGet(Field, getName());
    out = Zero();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_Zero_hpp_
