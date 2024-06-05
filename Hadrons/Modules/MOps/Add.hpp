#ifndef Hadrons_MOps_Add_hpp_
#define Hadrons_MOps_Add_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Add                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MOps)

class AddPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(AddPar,
                                    std::string, left,
                                    double,      leftfactor,
                                    std::string, right,
                                    double,      rightfactor);
};

template <typename Field>
class TAdd: public Module<AddPar>
{
public:
    // constructor
    TAdd(const std::string name);
    // destructor
    virtual ~TAdd(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(AddScalarField,     TAdd<FIMPL::Field          >, MOps);
MODULE_REGISTER_TMP(AddComplexField,    TAdd<FIMPL::ComplexField   >, MOps);
MODULE_REGISTER_TMP(AddFermionField,    TAdd<FIMPL::FermionField   >, MOps);
MODULE_REGISTER_TMP(AddPropagatorField, TAdd<FIMPL::PropagatorField>, MOps);

/******************************************************************************
 *                 TAdd implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TAdd<Field>::TAdd(const std::string name)
: Module<AddPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TAdd<Field>::getInput(void)
{
    std::vector<std::string> in {par().left, par().right};
    
    return in;
}

template <typename Field>
std::vector<std::string> TAdd<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TAdd<Field>::setup(void)
{
    if (!envHasType(Field, par().left))
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().left 
                          + "' has an unexpected type '"
                          + env().getObjectType(par().left)
                          + "' (expected " + typeName<Field>()
                          + ")", env().getObjectAddress(par().left))
    }
    if (!envHasType(Field, par().right))
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().right 
                          + "' has an unexpected type '"
                          + env().getObjectType(par().right)
                          + "' (expected " + typeName<Field>()
                          + ")", env().getObjectAddress(par().right))
    }

    unsigned int Ls_left  = env().getObjectLs(par().left );
    unsigned int Ls_right = env().getObjectLs(par().right);
    if (Ls_left != Ls_right)
    {
        std::string sError{ "Ls mismatch: left Ls="};
        sError.append( std::to_string( Ls_left ) );
        sError.append( ", right Ls=" );
        sError.append( std::to_string( Ls_right ) );
        HADRONS_ERROR(Size, sError);
    }
    envCreate(Field, getName(), Ls_left, envGetGrid(Field));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TAdd<Field>::execute(void)
{
    const Field& left  = envGet(Field, par().left);
    const Field& right = envGet(Field, par().right);
    Field&       out   = envGet(Field, getName());

    out = par().leftfactor * left + par().rightfactor * right;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MOps_Add_hpp_
