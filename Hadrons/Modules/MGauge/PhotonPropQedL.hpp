#ifndef Hadrons_MGauge_PhotonPropQedL_hpp_
#define Hadrons_MGauge_PhotonPropQedL_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         PhotonPropQedL                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class PhotonPropQedLPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PhotonPropQedLPar,
                                    std::string, improvement);
};

template <typename VType>
class TPhotonPropQedL: public Module<PhotonPropQedLPar>
{
public:
    typedef TEmFieldGenerator<VType>    EmGen;
    typedef typename EmGen::GaugeField  GaugeField;
    typedef typename EmGen::ScalarField PhotonProp;
public:
    // constructor
    TPhotonPropQedL(const std::string name);
    // destructor
    virtual ~TPhotonPropQedL(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(PhotonPropQedL, TPhotonPropQedL<vComplex>, MGauge);

/******************************************************************************
 *                 TPhotonPropQedL implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename VType>
TPhotonPropQedL<VType>::TPhotonPropQedL(const std::string name)
: Module<PhotonPropQedLPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename VType>
std::vector<std::string> TPhotonPropQedL<VType>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename VType>
std::vector<std::string> TPhotonPropQedL<VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename VType>
void TPhotonPropQedL<VType>::setup(void)
{
    envCreateLat(PhotonProp, getName());
    envTmp(EmGen, "gen", 1, envGetGrid(GaugeField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename VType>
void TPhotonPropQedL<VType>::execute(void)
{
    auto& a = envGet(PhotonProp, getName());
    std::vector<double> improvement = strToVec<double>(par().improvement);
    envGetTmp(EmGen, gen);
    
    LOG(Message) << "Generating momentum-space photon propagator (gauge: feynman)" << std::endl;
    if (improvement.size() > 0)
    {
        LOG(Message) << "Improvement coefficients " << improvement << std::endl;
    }
    gen.makeFeynmanPropQedL(a, improvement);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_PhotonPropQedL_hpp_
