#ifndef Hadrons_MUtilities_CShift_hpp_
#define Hadrons_MUtilities_CShift_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         CShift                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class CShiftPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CShiftPar,
                                    std::string, q,
                                    std::string, shift);
};

template <typename FImpl>
class TCShift: public Module<CShiftPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    // constructor
    TCShift(const std::string name);
    // destructor
    virtual ~TCShift(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CShift, TCShift<FIMPL>, MUtilities);

/******************************************************************************
 *                 TCShift implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCShift<FImpl>::TCShift(const std::string name)
: Module<CShiftPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCShift<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TCShift<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCShift<FImpl>::setup(void)
{
    envCreate(PropagatorField, getName(), 1, envGetGrid(PropagatorField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCShift<FImpl>::execute(void)
{
    // Set up shift variables
    Coordinate coord = strToVec<int>(par().shift);
    const PropagatorField& prop = envGet(PropagatorField, par().q);
    PropagatorField&        out = envGet(PropagatorField, getName());

    // Execute CShift
    out = prop;
    for (int mu=0;mu<coord.size();mu++) 
    {
        int shift = coord[mu];
        if (shift != 0)
            out = Cshift(out,mu,shift);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_CShift_hpp_
