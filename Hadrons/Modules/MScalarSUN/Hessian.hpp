#ifndef Hadrons_MScalarSUN_Hessian_hpp_
#define Hadrons_MScalarSUN_Hessian_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Hessian of a complex field                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class HessianPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HessianPar,
                                    std::string, op,
                                    std::string, dirs,
                                    DiffType,    type1,
                                    DiffType,    type2,
                                    std::string, output);
};

class HessianResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(HessianResult,
                                    DiffType,                          type1,
                                    DiffType,                          type2,
                                    std::vector<std::vector<Complex>>, value);
};

template <typename SImpl>
class THessian: public Module<HessianPar>
{
public:
    typedef typename SImpl::Field        Field;
    typedef typename SImpl::ComplexField ComplexField;
public:
    // constructor
    THessian(const std::string name);
    // destructor
    virtual ~THessian(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Hessian, THessian<SIMPL>, MScalarSUN);
MODULE_REGISTER_TMP(HessianSU2, THessian<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(HessianSU3, THessian<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(HessianSU4, THessian<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(HessianSU5, THessian<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(HessianSU6, THessian<ScalarNxNAdjImplR<6>>, MScalarSUN);

/******************************************************************************
 *                        THessian implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
THessian<SImpl>::THessian(const std::string name)
: Module<HessianPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> THessian<SImpl>::getInput(void)
{
    std::vector<std::string> in = {par().op};
    
    return in;
}

template <typename SImpl>
std::vector<std::string> THessian<SImpl>::getOutput(void)
{
    std::vector<std::string>  out;
    std::vector<unsigned int> dirs = strToVec<unsigned int>(par().dirs);

    for (auto &mu: dirs)
    for (auto &nu: dirs)
    {
        out.push_back(varName(getName(), mu, nu));
    }
    out.push_back(getName());

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void THessian<SImpl>::setup(void)
{
    std::vector<unsigned int> dirs = strToVec<unsigned int>(par().dirs);

    for (auto &mu: dirs)
    for (auto &nu: dirs)
    {
        envCreateLat(ComplexField, varName(getName(), mu, nu));
    }
    envTmpLat(ComplexField, "buf");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void THessian<SImpl>::execute(void)
{
    LOG(Message) << "Computing the " << par().type1 << "/" <<  par().type2 
        << " Hessian of '" << par().op << "'" << std::endl;
    
    std::vector<unsigned int> dirs = strToVec<unsigned int>(par().dirs);
    const unsigned int ndirs = dirs.size();
    HessianResult      result;
    auto               &op = envGet(ComplexField, par().op);

    if (!par().output.empty())
    {
        result.type1 = par().type1;
        result.type2 = par().type2;
        result.value.resize(ndirs, std::vector<Complex>(ndirs));
    }
    for (auto &nu: dirs)
    {
        envGetTmp(ComplexField, buf);

        dmu(buf, op, nu, par().type2);
        for (auto &mu: dirs)
        {
            auto &der = envGet(ComplexField, varName(getName(), mu, nu));
            
            dmu(der, buf, mu, par().type1);
            if (!par().output.empty())
            {
                result.value[mu][nu] = TensorRemove(sum(der));
            }
        }
    }

    saveResult(par().output, "hessian", result);
    envGet(HadronsSerializable, getName()) = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_Hessian_hpp_
