#ifndef Hadrons_MUtilities_CShift_hpp_
#define Hadrons_MUtilities_CShift_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

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

template <typename Field>
class TCShift: public Module<CShiftPar>
{
public:
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

    template<typename T>
    std::vector<T> getEnvVector(const std::string& value);
};

MODULE_REGISTER_TMP(CShiftGaugeField,      TCShift<GIMPL::GaugeField>,      MUtilities);
MODULE_REGISTER_TMP(CShiftPropagatorField, TCShift<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                 TCShift implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TCShift<Field>::TCShift(const std::string name)
: Module<CShiftPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TCShift<Field>::getInput(void)
{
    std::vector<std::string> in = {par().q};

    if ((!isVector<int>(par().shift)) && (!par().shift.empty()))
    {
        in.push_back(par().shift);
    }
    
    return in;
}

template <typename Field>
std::vector<std::string> TCShift<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TCShift<Field>::setup(void)
{
    envCreate(Field, getName(), 1, envGetGrid(Field));
}

template <typename Field>
template<typename T>
std::vector<T> TCShift<Field>::getEnvVector(const std::string& value)
{
    std::vector<T> res = strToVec<T>(value);
    if (res.size() <= 1)
    {
        if (envHasType(std::vector<T>, value))
        {
            LOG(Message) << "Trying to retrieve std::vector<T>" << std::endl;
            return envGet(std::vector<T>, value);
        }
        else if (envHasType(HadronsSerializable, value))
        {
            LOG(Message) << "Trying to retrieve HadronsSerializable" << std::endl;
            auto &t = envGet(HadronsSerializable, value);
            return t.template get<std::vector<T>>();
        }
        else
        {
            HADRONS_ERROR(Definition, "cannot interpret '" + value + "' as either a literal or environment-stored vector");
        }
    }
    LOG(Message) << "Trying to retrieve strToVec<T>" << std::endl;
    return res;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TCShift<Field>::execute(void)
{
    // Set up shift variables
    Coordinate             coord = getEnvVector<int>(par().shift);
    const Field& prop  = envGet(Field, par().q);
    Field&       out   = envGet(Field, getName());

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
