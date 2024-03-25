#ifndef Hadrons_MScalar_Loop_hpp_
#define Hadrons_MScalar_Loop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               Scalar loop                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class LoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoopPar,
                                    std::string, phi1,
                                    std::string, phi2,
                                    std::string, output);
};

template <typename SImpl>
class TLoop: public Module<LoopPar>
{
public:
    BASIC_TYPE_ALIASES(SImpl,);
public:
    // constructor
    TLoop(const std::string name);
    // destructor
    virtual ~TLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Loop, TLoop<SIMPL>, MScalar);

/******************************************************************************
 *                           TLoop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TLoop<SImpl>::TLoop(const std::string name)
: Module<LoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TLoop<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().phi1, par().phi2};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TLoop<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_field"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TLoop<SImpl>::setup(void)
{
    envCreateLat(ScalarField, getName() + "_field");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TLoop<SImpl>::execute(void)
{
    auto &phi1 = envGet(ScalarField, par().phi1);
    auto &phi2 = envGet(ScalarField, par().phi2);
    auto &g2   = envGet(ScalarField, getName() + "_field");

    g2 = phi1*phi2;

    std::vector<TComplex> buf;
    std::vector<Complex>  result;

    sliceSum(g2, buf, Tp);
    result.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result[t] = TensorRemove(buf[t]);
    }
    envGet(HadronsSerializable, getName()) = result;
    saveResult(par().output, "hvpLoop", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_Loop_hpp_
