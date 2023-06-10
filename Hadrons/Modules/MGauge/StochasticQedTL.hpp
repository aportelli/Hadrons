#ifndef Hadrons_MGauge_StochasticQedTL_hpp_
#define Hadrons_MGauge_StochasticQedTL_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Stochastic  QED_L/QED_TL                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class StochasticQedSubZmPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochasticQedSubZmPar,
                                    QedGauge, gauge);
};

template <typename VType>
class TStochasticQedTL: public Module<StochasticQedSubZmPar>
{
public:
    typedef TEmFieldGenerator<VType>    EmGen;
    typedef typename EmGen::GaugeField  GaugeField;
    typedef typename EmGen::ScalarField ScalarField;
public:
    // constructor
    TStochasticQedTL(const std::string name);
    // destructor
    virtual ~TStochasticQedTL(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool weightDone_;
};

MODULE_REGISTER_TMP(StochasticQedTL, ARG(TStochasticQedTL<vComplex>), MGauge);

/******************************************************************************
 *                 TStochasticQedTL implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename VType>
TStochasticQedTL<VType>::TStochasticQedTL(const std::string name)
: Module<StochasticQedSubZmPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename VType>
std::vector<std::string> TStochasticQedTL<VType>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename VType>
std::vector<std::string> TStochasticQedTL<VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename VType>
void TStochasticQedTL<VType>::setup(void)
{
    weightDone_ = env().hasCreatedObject("_" + getName() + "_weight");
    envCacheLat(ScalarField, "_" + getName() + "_weight");
    envCreateLat(GaugeField, getName());
    envTmp(EmGen, "gen", 1, envGetGrid(GaugeField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename VType>
void TStochasticQedTL<VType>::execute(void)
{
    auto &a = envGet(GaugeField, getName());
    auto &w = envGet(ScalarField, "_" + getName() + "_weight");
    envGetTmp(EmGen, gen);
    if (!weightDone_)
    {
        LOG(Message) << "Caching stochastic QED_TL EM potential weights" << std::endl;
        gen.makeWeightsQedTL(w);
    }
    LOG(Message) << "Generating stochastic EM potential (gauge: " << par().gauge << ")" << std::endl;
    auto tr = gen.getGaugeTranform(par().gauge);
    gen(a, rng4d(), w, tr);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_StochasticQedTL_hpp_
