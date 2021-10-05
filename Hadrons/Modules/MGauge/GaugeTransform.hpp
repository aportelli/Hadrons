#ifndef Hadrons_MGauge_GaugeTransform_hpp_
#define Hadrons_MGauge_GaugeTransform_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Gauge transformation on gauge field                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class GaugeTransformPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeTransformPar,
                                    std::string, gauge,
                                    std::string, transform);
};

template <typename GImpl>
class TGaugeTransform: public Module<GaugeTransformPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TGaugeTransform(const std::string name);
    // destructor
    virtual ~TGaugeTransform(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(GaugeTransform, TGaugeTransform<GIMPL>, MGauge);

/******************************************************************************
 *                 TGaugeTransform implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TGaugeTransform<GImpl>::TGaugeTransform(const std::string name)
: Module<GaugeTransformPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TGaugeTransform<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge, par().transform};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TGaugeTransform<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TGaugeTransform<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(GaugeLinkField, "adjg");
    envTmpLat(GaugeLinkField, "Umu");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TGaugeTransform<GImpl>::execute(void)
{
    LOG(Message) << "Gauge tranformation on gauge field '" << par().gauge
                 << "' using transformation matrix '" << par().transform << "'" << std::endl;

    auto &U   = envGet(GaugeField, par().gauge);
    auto &Utr = envGet(GaugeField, getName());
    auto &g   = envGet(GaugeLinkField, par().transform);
    envGetTmp(GaugeLinkField, adjg);
    envGetTmp(GaugeLinkField, Umu);

    adjg = adj(g);
    for(unsigned int mu = 0; mu < U.Grid()->Dimensions(); mu++)
    {
        Umu = PeekIndex<LorentzIndex>(U, mu);
        Umu = g*Umu*Cshift(adjg, mu, 1);
        PokeIndex<LorentzIndex>(Utr, Umu, mu);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_GaugeTransform_hpp_
