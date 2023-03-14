#ifndef Hadrons_MGauge_TopologicalCharge_hpp_
#define Hadrons_MGauge_TopologicalCharge_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         TopologicalCharge                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class TopologicalChargePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TopologicalChargePar,
                                    unsigned int, i);
};

template <typename FImpl>
class TTopologicalCharge: public Module<TopologicalChargePar>
{
public:
    // constructor
    TTopologicalCharge(const std::string name);
    // destructor
    virtual ~TTopologicalCharge(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TopologicalCharge, TTopologicalCharge<FIMPL>, MGauge);

/******************************************************************************
 *                 TTopologicalCharge implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTopologicalCharge<FImpl>::TTopologicalCharge(const std::string name)
: Module<TopologicalChargePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTopologicalCharge<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl>
std::vector<std::string> TTopologicalCharge<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTopologicalCharge<FImpl>::setup(void)
{

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTopologicalCharge<FImpl>::execute(void)
{

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_TopologicalCharge_hpp_
