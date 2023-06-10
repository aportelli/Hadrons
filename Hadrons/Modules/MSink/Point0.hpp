#ifndef Hadrons_MSink_Point0_hpp_
#define Hadrons_MSink_Point0_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                Point sink just picking the origin in space                 *
 *                                                                            *
 * e.g. sink(field(xvec, t)) = field(0, t)                                    * 
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

template <typename Field>
class TPoint0: public Module<NoPar>
{
public:
    typedef Field PropagatorField;
    typedef std::vector<typename PropagatorField::scalar_object> SlicedPropagator;
    typedef std::function<SlicedPropagator (const PropagatorField &)> SinkFn;

public:
    // constructor
    TPoint0(const std::string name);
    // destructor
    virtual ~TPoint0(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ScalarPoint0, TPoint0<ScalarImplCR::Field>, MSink);

/******************************************************************************
 *                            TPoint0 implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TPoint0<Field>::TPoint0(const std::string name)
: Module<NoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TPoint0<Field>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field>
std::vector<std::string> TPoint0<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TPoint0<Field>::setup(void)
{
    envCreate(SinkFn, getName(), 1, nullptr);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TPoint0<Field>::execute(void)
{
    LOG(Message) << "Setting up origin point sink function" << std::endl;
    auto sink = [this](const PropagatorField &field)
    {
        Coordinate origin;
        unsigned int nd = env().getNd();
        unsigned int nt = env().getDim(nd - 1);
        SlicedPropagator res(nt);

        for (unsigned int i = 0; i < nd - 1; ++i)
        {
            origin[i] = 0;
        }
        for (unsigned int t = 0; t < nt; ++t)
        {
            origin[nd - 1] = t;
            res[t] = peekSite(field, origin);
        }
        
        return res;
    };
    envGet(SinkFn, getName()) = sink;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_Point0_hpp_
