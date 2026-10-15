#ifndef Hadrons_MContraction_QEDSpecs_hpp_
#define Hadrons_MContraction_QEDSpecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         QEDSpecs                                 *
 ******************************************************************************/

/*
 Universal disconnected "Specs" loop subdiagram
             ___              ___
            /   \            /   \
           |     |~~~~~~~~~~|     |
            \___/   photon   \___/
              q               loop

 i Tr[q * Gamma_{mu}](x, x) * G^{mu,nu}(x, y) * i Tr[loop * Gamma_{nu}](y, y) 
 \__________________________________________/
                    Tadpole
*/


BEGIN_MODULE_NAMESPACE(MContraction)

class QEDSpecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(QEDSpecsPar,
                                    std::string, tadpole,
                                    std::string, loop,
                                    std::string, output)
};


template <typename FImpl, typename VType>
class TQEDSpecs: public Module<QEDSpecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);

    typedef TEmFieldGenerator<VType> EmGen;
    typedef typename EmGen::GaugeField  EmField;
    typedef typename EmGen::ScalarField PhotonProp;

    // constructor
    TQEDSpecs(const std::string name);
    // destructor
    virtual ~TQEDSpecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void) override;
    virtual std::vector<std::string> getOutput(void) override;
    virtual std::vector<std::string> getOutputFiles(void) override;
    // setup
    virtual void setup(void) override;
    // execution
    virtual void execute(void) override;
};

MODULE_REGISTER_TMP(QEDSpecs, ARG(TQEDSpecs<FIMPL, vComplex>), MContraction);

/******************************************************************************
 *                 TQEDSpecs implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
TQEDSpecs<FImpl, VType>::TQEDSpecs(const std::string name)
: Module<QEDSpecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename VType>
std::vector<std::string> TQEDSpecs<FImpl, VType>::getInput(void)
{
    std::vector<std::string> in = {par().tadpole, par().loop};
    
    return in;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDSpecs<FImpl, VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename VType>
std::vector<std::string> TQEDSpecs<FImpl, VType>::getOutputFiles(void)
{
    std::vector<std::string> output = {};
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDSpecs<FImpl, VType>::setup(void)
{
    envTmp(LatticeComplex, "tadpole_mu", 1, envGetGrid(ComplexField));
    envTmp(LatticeComplex, "looptrace",  1, envGetGrid(ComplexField));
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename VType>
void TQEDSpecs<FImpl, VType>::execute(void)
{
    Gamma::Algebra Gmu[] = 
    {
        (Gamma::Algebra::GammaX),
        (Gamma::Algebra::GammaY),
        (Gamma::Algebra::GammaZ),
        (Gamma::Algebra::GammaT),
    };

    // Fetch env variables
    LOG(Message) << "Starting Specs contraction using tadpole field '" << par().tadpole << "'" << std::endl;
    const EmField&         tadpole = envGet(EmField,         par().tadpole);
    const PropagatorField& loop    = envGet(PropagatorField, par().loop);
    envGetTmp(LatticeComplex, tadpole_mu);
    envGetTmp(LatticeComplex, looptrace);

    ComplexD result = 0;             // Output variable.

    // Perform the Lorentz index contraction between the
    // tadpole field and the EM current insertion on the loop.
    for (int mu=0; mu<env().getNd(); mu++)
    {
      looptrace  = trace(loop*Gamma(Gmu[mu]));
      tadpole_mu = peekLorentz(tadpole, mu);
      result += TensorRemove(sum(tadpole_mu * looptrace));
    }
    // The tadpole is assumed to already have the factor i from the
    // EM current insertion baked in. This i comes from the current
    // insertion on `loop`.
    result *= ComplexD(0.0,1.0);

    LOG(Message) << "specs: " << result << std::endl;

    saveResult(par().output, "specs", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_QEDSpecs_hpp_
