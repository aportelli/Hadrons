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
                                    std::string, field,
                                    std::string, output);
};

class PowerSpectrumResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PowerSpectrumResult,
                                    unsigned int, basisSize,
                                    unsigned int, vectorSize,
                                    std::vector<double>, eval,
                                    std::vector<std::vector<double>>, spectrum);
};

template <typename FImpl>
class TPowerSpectrum: public Module<PowerSpectrumPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
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

MODULE_REGISTER_TMP(FermionPowerSpectrum, TPowerSpectrum<FIMPL>, MCovariantFourier);

/******************************************************************************
 *                 TPowerSpectrum implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPowerSpectrum<FImpl>::TPowerSpectrum(const std::string name)
: Module<PowerSpectrumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPowerSpectrum<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().basis, par().field};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPowerSpectrum<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPowerSpectrum<FImpl>::setup(void)
{
    unsigned int Ls = env().getObjectLs(par().field);
    auto &field = envGet(BaseEigenPack<FermionField>, par().field);
    
    if (Ls > 1)
    {
        envTmp(ColourVectorField, "vec5", Ls, field.evec[0].Grid());
    }
    envTmp(ColourVectorField, "tmp", Ls, field.evec[0].Grid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPowerSpectrum<FImpl>::execute(void)
{
    auto                &basis = envGet(BaseEigenPack<ColourVectorField>, par().basis);
    auto                &field = envGet(BaseEigenPack<FermionField>, par().field);
    unsigned int        Ls     = env().getObjectLs(par().field);
    Real                coeff;
    PowerSpectrumResult res;
    envGetTmp(ColourVectorField, tmp);

    res.basisSize  = basis.evec.size();
    res.vectorSize = field.evec.size();
    res.eval.resize(res.basisSize);
    res.spectrum.resize(res.vectorSize, std::vector<double>(res.basisSize));
    for (unsigned int i = 0; i < field.evec.size(); ++i)
    {
        LOG(Message) << "vector " << i << " spectrum calculation" << std::endl;
        for (unsigned int j = 0; j < basis.evec.size(); ++j)
        {
            coeff = 0.;
            if (Ls == 1)
            {
                for (unsigned int s = 0; s < Ns; ++s)
                {
                    tmp = peekSpin(field.evec[i], s);
                    conformable(tmp, basis.evec[j]);
                    coeff += std::norm(innerProduct(tmp, basis.evec[j]));
                }
            }
            else
            {
                envGetTmp(ColourVectorField, vec5);

                for (unsigned int s = 0; s < Ls; ++s)
                {
                    InsertSlice(basis.evec[j], vec5, s, 0);
                }
                for (unsigned int s = 0; s < Ns; ++s)
                {
                    tmp = peekSpin(field.evec[i], s);
                    conformable(tmp, vec5);
                    coeff += std::norm(innerProduct(tmp, vec5));
                }
            }
            if (i == 0)
            {
                res.eval[j] = basis.eval[j];
            }
            res.spectrum[i][j] = coeff;
        }
    }
    saveResult(par().output, "spectrum", res);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MCovariantFourier_PowerSpectrum_hpp_
