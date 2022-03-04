#ifndef Hadrons_MWavelet_Spectrum_hpp_
#define Hadrons_MWavelet_Spectrum_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DWT.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Spectrum                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MWavelet)

class SpectrumPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SpectrumPar,
                                    std::string, gauge,
                                    std::string, field,
                                    std::string, output,
                                    std::string, filter,
                                    std::string, compressedPack,
                                    std::string, op,
                                    unsigned int, compCut,
                                    unsigned int, maxLevel);
};

class SpectrumResult: Serializable
{
public:
    typedef Eigen::Tensor<Real, 3> Tensor3;
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SpectrumResult,
                                    std::string, filter,
                                    std::vector<double>, norm,
                                    Tensor3, spectrum);
};

template <typename Field, typename GImpl>
class TSpectrum: public Module<SpectrumPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
    typedef LatticeWt<Field>           Wt;
    typedef CovariantDwt<Field, GImpl> CDwt;
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

MODULE_REGISTER_TMP(Spectrum, ARG(TSpectrum<FIMPL::FermionField, GIMPL>), MWavelet);
MODULE_REGISTER_TMP(ColourVectorSpectrum, ARG(TSpectrum<ColourVectorField<FIMPL>, GIMPL>), MWavelet);

/******************************************************************************
 *                 TSpectrum implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
TSpectrum<Field, GImpl>::TSpectrum(const std::string name)
: Module<SpectrumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field, typename GImpl>
std::vector<std::string> TSpectrum<Field, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().field, par().gauge};
    
    if (!par().compressedPack.empty())
    {
        in.push_back(par().op);
    }
    

    return in;
}

template <typename Field, typename GImpl>
std::vector<std::string> TSpectrum<Field, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    if (!par().compressedPack.empty())
    {
        out.push_back(par().compressedPack);
    }
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TSpectrum<Field, GImpl>::setup(void)
{
    envCache(std::vector<std::shared_ptr<GridBase>>, "gridCache", 1, 0);
    envCache(std::vector<GridBase *>, "gridPt", 1, 0);
    auto &gridCache = envGet(std::vector<std::shared_ptr<GridBase>>, "gridCache");
    auto &gridPt    = envGet(std::vector<GridBase *>, "gridPt");
    auto &field     = envGet(BaseEigenPack<Field>, par().field);
    
    if (gridCache.empty())
    {
        gridPt.push_back(getGrid<Field>());
        for (unsigned int l = 1; l <= par().maxLevel; ++l)
        {
            gridCache.emplace_back(Dwt<Field>::coarsenGrid(gridPt[l - 1]));
            gridPt.push_back(gridCache.back().get());
        }
    }
    envTmp(Wt, "wt", 1, gridPt, par().maxLevel);
    envGetTmp(Wt, wt);
    auto &U  = envGet(GaugeField, par().gauge);
    envTmp(CDwt, "dwt", 1, wt.getGrid(), *DwtFilters::fromName.at(par().filter), 
           U, par().maxLevel);
    envTmpLat(Field, "rec");
    envTmpLat(Field, "diff");
    if (!par().compressedPack.empty())
    {
        envCreateDerived(BaseEigenPack<Field>, EigenPack<Field>, 
            par().compressedPack, 1, field.evec.size(), getGrid<Field>());
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TSpectrum<Field, GImpl>::execute(void)
{
    SpectrumResult     result;
    auto               &field = envGet(BaseEigenPack<Field>, par().field);
    const unsigned int nd     = field.evec[0].Grid()->Nd();
    unsigned int       nComp  = 1.;
    
    envGetTmp(Wt, wt);
    envGetTmp(CDwt, dwt);
    envGetTmp(Field, rec);
    envGetTmp(Field, diff);
    for (unsigned int mu = 0; mu < nd; ++mu)
    {
        nComp *= 2;
    }
    result.spectrum.resize(field.evec.size(), par().maxLevel, nComp);
    result.filter = par().filter;
    if (!par().compressedPack.empty())
    {
        auto &pack = envGet(BaseEigenPack<Field>, par().compressedPack);

        pack.eval = field.eval;
    }
    for (unsigned int i = 0; i < field.evec.size(); ++i)
    {
        Real n2 = norm2(field.evec[i]);

        LOG(Message) << "----- vector " << i << std::endl;
        result.norm.push_back(n2);
        dwt.forward(wt, field.evec[i]);
        for (unsigned int l = 0; l < par().maxLevel; ++l)
        for (unsigned int j = 0; j < wt(l).size(); ++j)
        {
            result.spectrum(i, l, j) = norm2(wt(l, j));
            LOG(Message) << "Component (" << l << ", " << j << "): " << result.spectrum(i, l, j) << std::endl;
        }
        if (!par().compressedPack.empty())
        {
            auto &pack = envGet(BaseEigenPack<Field>, par().compressedPack);
            auto &op   = envGet(LinearOperatorBase<Field>, par().op);

            for (unsigned int j = 1; j < wt(0).size(); ++j)
            {
                size_t numHigh = 0;

                for(size_t k = 0; k < CHAR_BIT * sizeof(j); ++k)
                {
                if ((j & (1 << k)) != 0)
                    ++numHigh;
                }

                LOG(Message) << j << " " << numHigh << std::endl;
                if (numHigh >= par().compCut)
                {
                    wt(0, j) = Zero();
                }
            }
            dwt.backward(pack.evec[i], wt);

            Field diff(getGrid<Field>()), tmp(getGrid<Field>());

            diff = pack.evec[i] - field.evec[i];
            LOG(Message) << "norm2 diff= " << norm2(diff) << std::endl;
            op.HermOp(field.evec[i], tmp);
            diff = tmp - field.eval[i]*field.evec[i];
            LOG(Message) << "norm2 EV check= " << norm2(diff) << std::endl;
            op.HermOp(pack.evec[i], tmp);
            diff = tmp - pack.eval[i]*pack.evec[i];
            LOG(Message) << "norm2 EV diff= " << norm2(diff) << std::endl;
            for (unsigned int j = 1; j < field.evec.size(); ++j)
            {
                LOG(Message) << "(ci, " << j << ")= " << innerProduct(pack.evec[i], field.evec[j]) << std::endl;
            }
        }
    }
    saveResult(par().output, "spectrum", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MWavelet_Spectrum_hpp_
