#ifndef Hadrons_MGauge_ExactQedL_hpp_
#define Hadrons_MGauge_ExactQedL_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ExactQedL                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class ExactQedLPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactQedLPar,
                                    QedGauge, gauge,
                                    std::string, improvement);
};

template <typename VType>
class TExactQedL: public Module<ExactQedLPar>
{
public:
    typedef TEmFieldGenerator<VType>    EmGen;
    typedef typename EmGen::GaugeField  GaugeField;
    typedef typename EmGen::ScalarField ScalarField;
public:
    // constructor
    TExactQedL(const std::string name);
    // destructor
    virtual ~TExactQedL(void) {};
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

MODULE_REGISTER_TMP(ExactQedL, TExactQedL<vComplex>, MGauge);

/******************************************************************************
 *                 TExactQedL implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename VType>
TExactQedL<VType>::TExactQedL(const std::string name)
: Module<ExactQedLPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename VType>
std::vector<std::string> TExactQedL<VType>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename VType>
std::vector<std::string> TExactQedL<VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename VType>
void TExactQedL<VType>::setup(void)
{
    weightDone_ = env().hasCreatedObject("_" + getName() + "_weight");
    envCacheLat(ScalarField, "_" + getName() + "_weight");
    envCreateLat(GaugeField, getName());
    envTmp(EmGen, "gen", 1, envGetGrid(GaugeField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename VType>
void TExactQedL<VType>::execute(void)
{
    auto &a = envGet(GaugeField, getName());
    auto &w = envGet(ScalarField, "_" + getName() + "_weight");
    std::vector<double> improvement = strToVec<double>(par().improvement);
    envGetTmp(EmGen, gen);
    if (!weightDone_)
    {
        LOG(Message) << "Caching exact QED_L EM potential weights" << std::endl;
        if (improvement.size() > 0)
        {
            LOG(Message) << "Improvement coefficients " << improvement << std::endl;
        }
        gen.makeWeightsQedL(w, improvement);
    }
    LOG(Message) << "Generating exact EM potential (gauge: " << par().gauge << ")" << std::endl;
    auto tr = gen.getGaugeTranform(par().gauge);
    gen(a, w, tr);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_ExactQedL_hpp_
