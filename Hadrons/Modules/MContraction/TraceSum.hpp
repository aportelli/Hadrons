#ifndef Hadrons_MContraction_TraceSum_hpp_
#define Hadrons_MContraction_TraceSum_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>


BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         TraceSum                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class TraceSumPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TraceSumPar,
                                    std::string, q,
                                    Gamma::Algebra, gamma);
};

template <typename FImpl>
class TTraceSum: public Module<TraceSumPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);

    // constructor
    TTraceSum(const std::string name);
    // destructor
    virtual ~TTraceSum(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TraceSum, TTraceSum<FIMPL>, MContraction);

/******************************************************************************
 *                 TTraceSum implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTraceSum<FImpl>::TTraceSum(const std::string name)
: Module<TraceSumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTraceSum<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTraceSum<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_latticeSum"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTraceSum<FImpl>::setup(void)
{
    envCreateLat(LatticeComplex, getName());
    envCreate(HadronsSerializable, getName() + "_latticeSum", 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTraceSum<FImpl>::execute(void)
{   
    const auto& q = envGet(PropagatorField, par().q);
    auto& out     = envGet(LatticeComplex, getName());
    auto& out_sum = envGet(HadronsSerializable, getName() + "_latticeSum");

    Gamma gamma(par().gamma);
    out     = trace(q*gamma);
    out_sum = sum(out);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_TraceSum_hpp_
