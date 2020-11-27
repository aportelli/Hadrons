#ifndef Hadrons_MScalarSUN_Reweight_hpp_
#define Hadrons_MScalarSUN_Reweight_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Reweight                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class ReweightPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ReweightPar,
                                    RealD      , m2,
                                    RealD      , lambda,
                                    RealD      , g,
                                    std::string, field,
                                    std::string,  output);
};

class ReweightResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ReweightResult,
                                    Real, value);
};

template <typename FImpl>
class TReweight: public Module<ReweightPar>
{
public:
    typedef typename FImpl::Field        Field;
    typedef typename FImpl::ComplexField ComplexField;
public:
    // constructor
    TReweight(const std::string name);
    // destructor
    virtual ~TReweight(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ReweightSU2, TReweight<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(ReweightSU3, TReweight<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(ReweightSU4, TReweight<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(ReweightSU5, TReweight<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(ReweightSU6, TReweight<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                 TReweight implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TReweight<FImpl>::TReweight(const std::string name)
: Module<ReweightPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TReweight<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    return in;
}

template <typename FImpl>
std::vector<std::string> TReweight<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TReweight<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TReweight<FImpl>::execute(void)
{
    LOG(Message) << "Computing reweight factor..." << std::endl;
    LOG(Message) << "g = " << par().g << std::endl;
    LOG(Message) << "m2 = " << par().m2 << std::endl;
    LOG(Message) << "lambda = " << par().lambda << std::endl;
    
    ReweightResult result;
    auto           &phi = envGet(Field, par().field);
    const unsigned int N = FImpl::Group::Dimension, nd = env().getNd();
    
    ComplexField Trphi = trace(phi);
    ComplexField Trphi2 = trace(phi*phi);
    ComplexField Trphi3 = trace(phi*phi*phi);
    ComplexField Tr2dphi = 0.*trace(phi);
    
    for (int mu=0;mu<nd;++mu) Tr2dphi = Tr2dphi + trace(Cshift(phi,mu,1)-phi)*trace(Cshift(phi,mu,1)-phi);

    result.value = (1./par().g)*real(TensorRemove(sum( Tr2dphi + par().m2*Trphi*Trphi + par().lambda*(  (-3./pow(N,2))*Trphi*Trphi*Trphi*Trphi + (6./N)*Trphi*Trphi*Trphi2 - 4.*Trphi*Trphi3  )  )));


    LOG(Message) << "Trace(phi) sum = " << TensorRemove(sum(timesMinusI(Trphi))).real() << std::endl;
    LOG(Message) << "Trace^2(dphi) sum = " << -TensorRemove(sum(Tr2dphi)).real() << std::endl;
    LOG(Message) << "Trace(phi^2) sum = " << -TensorRemove(sum(Trphi2)).real() << std::endl;
    LOG(Message) << "Trace(phi^3) sum = " << -TensorRemove(sum(timesI(Trphi3))).real() << std::endl;
    LOG(Message) << "Reweight Factor = " << result.value << std::endl;
    saveResult(par().output, "reweight", result);
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Reweight_hpp_
