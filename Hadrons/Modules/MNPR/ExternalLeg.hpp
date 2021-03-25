#ifndef Hadrons_MNPR_ExternalLeg_hpp_
#define Hadrons_MNPR_ExternalLeg_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ExternalLeg                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class ExternalLegPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExternalLegPar,
                                    std::string, prop,
                                    std::string, momentum,
                                    bool, outgoing,
                                    std::string, output);
};

template <typename FImpl>
class TExternalLeg: public Module<ExternalLegPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,)
    // constructor
    TExternalLeg(const std::string name);
    // destructor
    virtual ~TExternalLeg(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExternalLeg, TExternalLeg<FIMPL>, MNPR);

/******************************************************************************
 *                 TExternalLeg implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TExternalLeg<FImpl>::TExternalLeg(const std::string name)
: Module<ExternalLegPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TExternalLeg<FImpl>::getInput(void)
{
    std::vector<std::string> in = { par().prop };

    return in;
}

template <typename FImpl>
std::vector<std::string> TExternalLeg<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExternalLeg<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for external leg '"
        << getName() << "'" << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TExternalLeg<FImpl>::execute(void)
{
    std::vector<int> latt_size(env().getGrid()->FullDimensions().toVector());

    PropagatorField prop_unphased = envGet(PropagatorField, par().prop);
    LatticeSpinColourMatrix prop_phased(env().getGrid());

    std::vector<Real> momentum = strToVec<Real>(par().momentum);

    LatticeComplex pdotx(env().getGrid());
    LatticeComplex coordinate(env().getGrid());

    RealD volume = 1.0;
    for (int mu = 0; mu < Nd; mu++) {
        volume *= latt_size[mu];
    }

    LOG(Message) << "Phasing propagator" << std::endl;

    pdotx = Zero();
    for (int mu = 0; mu < Nd; mu++) {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;
        pdotx += coordinate * momentum[mu];
    }
    Complex imag = Complex(0, -1.0);
    prop_phased = prop_unphased * exp(imag * pdotx);

    LOG(Message) << "Done phasing propagators" << std::endl;

    LOG(Message) << "Comptuing average" << std::endl;
    SpinColourMatrix average = (1.0 / volume) * sum(prop_phased);

    if (par().outgoing) {
        // Reversing quark flow for outgoing propagators
        Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
        average = g5 * adj(average) * g5;
    }
    LOG(Message) << "Done comptuing average" << std::endl;

    LOG(Message) << "Saving output to '" << par().output << "'" << std::endl;
    saveResult(par().output, "external_leg", average);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNPR_ExternalLeg_hpp_
