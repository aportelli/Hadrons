#ifndef Hadrons_MCovariantFourier_PowerSpectrum_hpp_
#define Hadrons_MCovariantFourier_PowerSpectrum_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         PowerSpectrum                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MCovariantFourier)

class PowerSpectrumPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PowerSpectrumPar,
                                    std::string, basis,
                                    std::string, field);
};

template <typename Field>
class TPowerSpectrum: public Module<PowerSpectrumPar>
{
public:
    // constructor
    TPowerSpectrum(const std::string name);
    // destructor
    virtual ~TPowerSpectrum(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionPowerSpectrum, 
                    TPowerSpectrum<FIMPL::FermionField>, MCovariantFourier);

/******************************************************************************
 *                 TPowerSpectrum implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TPowerSpectrum<Field>::TPowerSpectrum(const std::string name)
: Module<PowerSpectrumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TPowerSpectrum<Field>::getInput(void)
{
    std::vector<std::string> in = {par().basis, par().field};
    
    return in;
}

template <typename Field>
std::vector<std::string> TPowerSpectrum<Field>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TPowerSpectrum<Field>::setup(void)
{
    unsigned int Ls = env().getObjectLs(par().field);
    
    if (Ls > 1)
    {
        auto &field = envGet(BaseEigenPack<Field>, par().field);
        envTmp(Field, "vec5", Ls, field.evec[0].Grid());
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TPowerSpectrum<Field>::execute(void)
{
    auto    &epack = envGet(BaseEigenPack<Field>, par().basis);
    auto    &field = envGet(BaseEigenPack<Field>, par().field);
    unsigned int Ls = env().getObjectLs(par().field);
    Complex coeff;

    for (unsigned int i = 0; i < field.evec.size(); ++i)
    {
        for (unsigned int j = 0; j < epack.evec.size(); ++j)
        {
            envGetTmp(Field, vec5);

            for (unsigned int s = 0; s < Ls; ++s)
            {
                InsertSlice(epack.evec[j], vec5, s, 0);
            }
            conformable(vec5, field.evec[i]);
            coeff = innerProduct(vec5, field.evec[i]);
            LOG(Message) << "vector " << i << " mode " << j << " -- lambda= " 
                         << epack.eval[j] << " -- coeff= " << norm(coeff) << std::endl;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MCovariantFourier_PowerSpectrum_hpp_
