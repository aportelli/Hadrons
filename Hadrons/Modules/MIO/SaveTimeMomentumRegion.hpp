#ifndef Hadrons_MIO_SaveTimeMomentumRegion_hpp_
#define Hadrons_MIO_SaveTimeMomentumRegion_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SaveTimeMomentumRegion                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class SaveTimeMomentumRegionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveTimeMomentumRegionPar,
                                    std::string, tmomField,
                                    unsigned int, cut,
                                    std::string, output);
};

template <typename Field>
class TSaveTimeMomentumRegion: public Module<SaveTimeMomentumRegionPar>
{
public:
    typedef typename Field::scalar_object Site;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<std::vector<int>>, momentum,
                                        std::vector<std::vector<Site>>, corr);
    };
    class ResultVector: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ResultVector,
                                        std::vector<std::vector<int>>, momentum,
                                        std::vector<std::vector<std::vector<Site>>>, corr);
    };
public:
    // constructor
    TSaveTimeMomentumRegion(const std::string name);
    // destructor
    virtual ~TSaveTimeMomentumRegion(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    // row-major momentum indexing
    static void indexToMom(std::vector<int> &mom, const unsigned int i,
                           const unsigned int cut, const unsigned int nd);
};

MODULE_REGISTER_TMP(SaveScalarTimeMomentumRegion, TSaveTimeMomentumRegion<SIMPL::ComplexField>, MIO);

/******************************************************************************
 *                 TSaveTimeMomentumRegion implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TSaveTimeMomentumRegion<Field>::TSaveTimeMomentumRegion(const std::string name)
: Module<SaveTimeMomentumRegionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TSaveTimeMomentumRegion<Field>::getInput(void)
{
    std::vector<std::string> in = {par().tmomField};
    
    return in;
}

template <typename Field>
std::vector<std::string> TSaveTimeMomentumRegion<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename Field>
std::vector<std::string> TSaveTimeMomentumRegion<Field>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
    {
        output.push_back(resultFilename(par().output));
    }
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TSaveTimeMomentumRegion<Field>::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TSaveTimeMomentumRegion<Field>::indexToMom(std::vector<int> &mom, 
                                                const unsigned int i,
                                                const unsigned int cut,
                                                const unsigned int nd)
{
    unsigned int buf, prod;
    const unsigned int size = 2*cut + 1;

    if (mom.size() != nd)
    {
        mom.resize(nd);
    }
    buf = i;
    prod = 1;
    for (int mu = nd - 1; mu >= 0; --mu)
    {
        mom[mu]  = (buf/prod)%size;
        buf     -= prod*mom[mu];
        mom[mu] -= nd;
        prod    *= size;
    }
}

template <typename Field>
void TSaveTimeMomentumRegion<Field>::execute(void)
{
    auto &field = envGet(Field, par().tmomField);
    const unsigned int nd = env().getNd(), nt = env().getDim(Tp), momSize = 2*par().cut + 1;
    unsigned int nMom = 1;
    std::vector<int> mom(nd - 1);
    Coordinate site(nd);
    Result result;

    for (unsigned int d = 0; d < nd - 1; ++d)
    {
        nMom *= momSize;
    }
    result.momentum.resize(nMom);
    result.corr.resize(nMom);
    for (unsigned int i = 0; i < nMom; ++i)
    {
        unsigned int j = 0;

        indexToMom(mom, i, par().cut, nd - 1);
        result.momentum[i] = mom;
        result.corr[i].resize(nt);
        for (unsigned int d = 0; d < nd; ++d)
        {
            if (d != Tp)
            {
                site[d] = mom[j];
                j++;
            }
        }
        for (unsigned int t = 0; t < nt; ++t)
        {
            site[Tp] = t;
            peekSite(result.corr[i][t], field, site);
        }
    }
    saveResult(par().output, "tmom_region", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveTimeMomentumRegion_hpp_
