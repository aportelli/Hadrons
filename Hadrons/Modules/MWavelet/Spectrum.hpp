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
                                    unsigned int, maxLevel);
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
    
    return in;
}

template <typename Field, typename GImpl>
std::vector<std::string> TSpectrum<Field, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
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
    
    if (gridCache.empty())
    {
        gridPt.push_back(getGrid<Field>());
        for (unsigned int l = 1; l <= par().maxLevel; ++l)
        {
            gridCache.emplace_back(Dwt<Field>::coarsenGrid(gridPt[l - 1]));
            gridPt.push_back(gridCache.back().get());
        }
    }
    envCreate(Wt, "wt", 1, gridPt, par().maxLevel);
    auto &wt = envGet(Wt, "wt");
    auto &U  = envGet(GaugeField, par().gauge);
    envCreate(CDwt, "dwt", 1, wt.getGrid(), 
              *DwtFilters::fromName.at(par().filter), U, par().maxLevel);
    envTmpLat(Field, "rec");
    envTmpLat(Field, "diff");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TSpectrum<Field, GImpl>::execute(void)
{
    auto &field = envGet(BaseEigenPack<Field>, par().field);
    auto &wt    = envGet(Wt, "wt");
    auto &dwt   = envGet(CDwt, "dwt");
    envGetTmp(Field, rec);
    envGetTmp(Field, diff);

    for (unsigned int i = 0; i < field.evec.size(); ++i)
    {
        LOG(Message) << "----- vector " << i << std::endl;
        dwt.forward(wt, field.evec[i]);
        for (unsigned int l = 0; l < par().maxLevel; ++l)
        {
            Real n2 = 1.;

            LOG(Message) << "Component (" << l << ", 0): " << norm2(wt(l, 0))/n2 << std::endl;
            for (unsigned int j = 1; j < wt(l).size(); ++j)
            {
                LOG(Message) << "Component (" << l << ", " << j << "): " << norm2(wt(l, j))/n2 << std::endl;
                wt(l, j) = Zero();
            }
            dwt.backward(rec, wt);
            diff = rec - field.evec[i];
            LOG(Message) << "Norm 2 relative diff= " << norm2(diff)/n2 << std::endl; 
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MWavelet_Spectrum_hpp_
