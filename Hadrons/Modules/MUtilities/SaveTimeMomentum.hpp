#ifndef Hadrons_MUtilities_SaveTimeMomentum_hpp_
#define Hadrons_MUtilities_SaveTimeMomentum_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SaveTimeMomentum                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class SaveTimeMomentumPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveTimeMomentumPar,
                                    std::string, tmomField,
                                    std::string, momentum,
                                    std::string, output);
};

template <typename Field>
class TSaveTimeMomentum: public Module<SaveTimeMomentumPar>
{
public:
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata, std::vector<int>, momentum);
    };
    typedef typename Field::scalar_object Site;
    typedef Correlator<Metadata, Site> Result;
public:
    // constructor
    TSaveTimeMomentum(const std::string name);
    // destructor
    virtual ~TSaveTimeMomentum(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveEmFieldTimeMomentum, TSaveTimeMomentum<TEmFieldGenerator<vComplex>::GaugeField>, MUtilities);

/******************************************************************************
 *                 TSaveTimeMomentum implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TSaveTimeMomentum<Field>::TSaveTimeMomentum(const std::string name)
: Module<SaveTimeMomentumPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TSaveTimeMomentum<Field>::getInput(void)
{
    std::vector<std::string> in = {par().tmomField};
    
    return in;
}

template <typename Field>
std::vector<std::string> TSaveTimeMomentum<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename Field>
std::vector<std::string> TSaveTimeMomentum<Field>::getOutputFiles(void)
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
void TSaveTimeMomentum<Field>::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TSaveTimeMomentum<Field>::execute(void)
{
    auto &field = envGet(Field, par().tmomField);
    unsigned int nd = env().getNd(), nt = env().getDim(Tp);
    Coordinate site(nd);
    std::vector<int> mom = strToVec<int>(par().momentum);
    Result result;

    if (mom.size() != env().getNd() - 1)
    {
        HADRONS_ERROR(Size, "momentum has " + std::to_string(mom.size())
                      + " components (must have " + std::to_string(env().getNd() - 1) + ")");
    }
    result.corr.resize(nt);
    unsigned int j = 0;
    Site f;
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
        peekSite(result.corr[t], field, site);
    }
    result.info.momentum = mom;
    saveResult(par().output, "meson", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_SaveTimeMomentum_hpp_
