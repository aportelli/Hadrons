#ifndef Hadrons_MGauge_StochasticQedSubZm_hpp_
#define Hadrons_MGauge_StochasticQedSubZm_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         StochasticQedSubZm                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class StochasticQedSubZmPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochasticQedSubZmPar,
                                    QedGauge, gauge);
};

enum ZmScheme { qedL = 0, qedTL = 1 };

template <typename VType, ZmScheme scheme>
class TStochasticQedSubZm: public Module<StochasticQedSubZmPar>
{
public:
    typedef TEmFieldGenerator<VType>    EmGen;
    typedef typename EmGen::GaugeField  GaugeField;
    typedef typename EmGen::ScalarField ScalarField;
public:
    // constructor
    TStochasticQedSubZm(const std::string name);
    // destructor
    virtual ~TStochasticQedSubZm(void) {};
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

MODULE_REGISTER_TMP(StochasticQedL, ARG(TStochasticQedSubZm<vComplex, ZmScheme::qedL>), MGauge);
MODULE_REGISTER_TMP(StochasticQedTL, ARG(TStochasticQedSubZm<vComplex, ZmScheme::qedTL>), MGauge);

/******************************************************************************
 *                 TStochasticQedSubZm implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename VType, ZmScheme scheme>
TStochasticQedSubZm<VType, scheme>::TStochasticQedSubZm(const std::string name)
: Module<StochasticQedSubZmPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename VType, ZmScheme scheme>
std::vector<std::string> TStochasticQedSubZm<VType, scheme>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename VType, ZmScheme scheme>
std::vector<std::string> TStochasticQedSubZm<VType, scheme>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename VType, ZmScheme scheme>
void TStochasticQedSubZm<VType, scheme>::setup(void)
{
    weightDone_ = env().hasCreatedObject("_" + getName() + "_weight");
    envCacheLat(ScalarField, "_" + getName() + "_weight");
    envCreateLat(GaugeField, getName());
    envTmp(EmGen, "gen", 1, envGetGrid(GaugeField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename VType, ZmScheme scheme>
void TStochasticQedSubZm<VType, scheme>::execute(void)
{
    auto &a = envGet(GaugeField, getName());
    auto &w = envGet(ScalarField, "_" + getName() + "_weight");
    envGetTmp(EmGen, gen);
    if (!weightDone_)
    {
        ;
        switch (scheme)
        {
        case ZmScheme::qedL:
            LOG(Message) << "Caching stochastic EM potential weights, zero-mode scheme: QED_L" << std::endl;
            gen.makeWeightsQedL(w);
            break;
        case ZmScheme::qedTL:
            LOG(Message) << "Caching stochastic EM potential weights, zero-mode scheme: QED_TL" << std::endl;
            gen.makeWeightsQedTL(w);
            break;
        default:
            HADRONS_ERROR(Definition, "invalid zero-mode scheme")
            break;
        }
    }
    LOG(Message) << "Generating stochastic EM potential (gauge: " << par().gauge << ")" << std::endl;
    auto tr = gen.getGaugeTranform(par().gauge);
    gen(a, rng4d(), w, tr);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_StochasticQedSubZm_hpp_
