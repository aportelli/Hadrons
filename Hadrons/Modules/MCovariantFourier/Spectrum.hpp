#ifndef Hadrons_MCovariantFourier_Spectrum_hpp_
#define Hadrons_MCovariantFourier_Spectrum_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Spectrum                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MCovariantFourier)

class SpectrumPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SpectrumPar,
                                    std::string, basis,
                                    std::string, field,
                                    std::string, output);
};

class SpectrumResult: Serializable
{
public:
    typedef Eigen::Tensor<Real, 2>    Tensor2;
    typedef Eigen::Tensor<Complex, 4> Tensor4;
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SpectrumResult,
                                    unsigned int, basisSize,
                                    unsigned int, vectorSize,
                                    std::vector<double>, eval,
                                    Tensor2, power,
                                    Tensor4, spectrum);
};

template <typename FImpl>
class TSpectrum: public Module<SpectrumPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSpectrum(const std::string name);
    // destructor
    virtual ~TSpectrum(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionSpectrum, TSpectrum<FIMPL>, MCovariantFourier);

/******************************************************************************
 *                 TSpectrum implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSpectrum<FImpl>::TSpectrum(const std::string name)
: Module<SpectrumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSpectrum<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().basis, par().field};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSpectrum<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSpectrum<FImpl>::setup(void)
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
void TSpectrum<FImpl>::execute(void)
{
    auto                &basis = envGet(BaseEigenPack<ColourVectorField>, par().basis);
    auto                &field = envGet(BaseEigenPack<FermionField>, par().field);
    unsigned int        Ls     = env().getObjectLs(par().field);
    SpectrumResult res;
    Real                energy, bandEnergy, n2;
    envGetTmp(ColourVectorField, tmp);

    res.basisSize  = basis.evec.size();
    res.vectorSize = field.evec.size();
    res.eval.resize(res.basisSize);
    res.power.resize(res.vectorSize, res.basisSize);
    res.power.setZero();
    res.spectrum.resize(res.vectorSize, res.basisSize, Ns, Ls);
    for (unsigned int i = 0; i < field.evec.size(); ++i)
    {
        LOG(Message) << "vector " << i << " spectrum calculation" << std::endl;
        energy     = norm2(field.evec[i]);
        bandEnergy = 0.;
        if (i == 0)
        {
            res.eval = basis.eval;
        }
        if (Ls == 1)
        {
            for (unsigned int s = 0; s < Ns; ++s)
            {
                tmp = peekSpin(field.evec[i], s);
                LOG(Message) << "spin= " << s << " / norm2= " << norm2(tmp) << std::endl;
                for (unsigned int j = 0; j < basis.evec.size(); ++j)
                {
                    conformable(tmp, basis.evec[j]);
                    res.spectrum(i, j, s, 0)  = innerProduct(tmp, basis.evec[j]);
                    n2                        = std::norm(res.spectrum(i, j, s, 0));
                    res.power(i, j)          += n2;
                    bandEnergy               += n2;
                }
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
                    LOG(Message) << "spin= " << s << " / s= " << t 
                                 << " / norm2= " << norm2(vec4) << std::endl;
                    for (unsigned int j = 0; j < basis.evec.size(); ++j)
                    {
                        conformable(vec4, basis.evec[j]);
                        res.spectrum(i, j, s, t)  = innerProduct(vec4, basis.evec[j]);
                        n2                        = std::norm(res.spectrum(i, j, s, t));
                        res.power(i, j)          += n2;
                        bandEnergy               += n2;
                    }
                }  
            }
        }
        LOG(Message) << "energy= " << energy << " / band energy= " << bandEnergy << " / lossyness= " << 1. - sqrt(bandEnergy/energy) << std::endl;
    }
    saveResult(par().output, "spectrum", res);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MCovariantFourier_Spectrum_hpp_
