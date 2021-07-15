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
    typedef Eigen::Tensor<Real, 2>    Tensor2;
    typedef Eigen::Tensor<Complex, 4> Tensor4;
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PowerSpectrumResult,
                                    unsigned int, basisSize,
                                    unsigned int, vectorSize,
                                    std::vector<double>, eval,
                                    Tensor2, power,
                                    Tensor4, spectrum);
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
        auto &basis = envGet(BaseEigenPack<ColourVectorField>, par().basis);
        envTmp(ColourVectorField, "vec4", 1, basis.evec[0].Grid());
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
    PowerSpectrumResult res;
    Real                energy, bandEnergy;
    envGetTmp(ColourVectorField, tmp);

    res.basisSize  = basis.evec.size();
    res.vectorSize = field.evec.size();
    res.eval.resize(res.basisSize);
    res.power.resize(res.vectorSize, res.basisSize);
    res.spectrum.resize(res.vectorSize, res.basisSize, Ns, Ls);
    for (unsigned int i = 0; i < field.evec.size(); ++i)
    {
        LOG(Message) << "vector " << i << " spectrum calculation" << std::endl;
        energy     = norm2(field.evec[i]);
        bandEnergy = 0.;
        for (unsigned int j = 0; j < basis.evec.size(); ++j)
        {
            res.power(i, j) = 0.;
            if (Ls == 1)
            {
                for (unsigned int s = 0; s < Ns; ++s)
                {
                    tmp = peekSpin(field.evec[i], s);
                    conformable(tmp, basis.evec[j]);
                    res.spectrum(i, j, s, 0)  = innerProduct(tmp, basis.evec[j]);
                    res.power(i, j)          += std::norm(res.spectrum(i, j, s, 0));
                    bandEnergy               += res.power(i, j) ;
                }
            }
            else
            {
                envGetTmp(ColourVectorField, vec4);
                for (unsigned int s = 0; s < Ns; ++s)
                {
                    tmp = peekSpin(field.evec[i], s);
                    for (unsigned int t = 0; t < Ls; ++t)
                    {
                        ExtractSlice(vec4, tmp, t, 0);
                        conformable(vec4, basis.evec[j]);
                        res.spectrum(i, j, s, t)  = innerProduct(vec4, basis.evec[j]);
                        res.power(i, j)          += std::norm(res.spectrum(i, j, s, t));
                        bandEnergy               += res.power(i, j) ;
                    }  
                }
            }
            if (i == 0)
            {
                res.eval[j] = basis.eval[j];
            }
        }
        LOG(Message) << "energy= " << energy << " / band energy= " << bandEnergy << " / lossyness= " << 1. - sqrt(bandEnergy/energy) << std::endl;
    }
    saveResult(par().output, "spectrum", res);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MCovariantFourier_PowerSpectrum_hpp_
