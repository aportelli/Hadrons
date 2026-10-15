#ifndef Hadrons_MUtilities_FFT_hpp_
#define Hadrons_MUtilities_FFT_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       FFT with dimension selection                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class FFTPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FFTPar,
                                    std::string, field,
                                    std::string, dimMask,
                                    bool, backward);
};

template <typename Field>
class TFFT: public Module<FFTPar>
{
public:
    // constructor
    TFFT(const std::string name);
    // destructor
    virtual ~TFFT(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void doFft(Field &out, const Field &in, const std::vector<int> &mask, 
               bool backward);
private:
    bool isVector_;
};

MODULE_REGISTER_TMP(ScalarFFT, TFFT<SIMPL::ComplexField>, MUtilities);
MODULE_REGISTER_TMP(EmFieldFFT, TFFT<TEmFieldGenerator<vComplex>::GaugeField>, MUtilities);

/******************************************************************************
 *                           TFFT implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TFFT<Field>::TFFT(const std::string name)
: Module<FFTPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TFFT<Field>::getInput(void)
{
    std::vector<std::string> in = {par().field};
    
    return in;
}

template <typename Field>
std::vector<std::string> TFFT<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TFFT<Field>::setup(void)
{
    isVector_ = envHasType(std::vector<Field>, par().field);
    if (isVector_)
    {
        auto &field = envGet(std::vector<Field>, par().field);
        if (field.size() == 0)
        {
            HADRONS_ERROR(Size, "input field vector is empty");
        }
        envCreate(std::vector<Field>, getName(), 1, field.size(), field[0].Grid());
    }
    else
    {
        envCreateLat(Field, getName());
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
inline void TFFT<Field>::doFft(Field &out, const Field &in, 
                               const std::vector<int> &mask, bool backward)
{
    GridBase *g = in.Grid();
    FFT fft(dynamic_cast<GridCartesian *>(g));
    Coordinate maskc(mask);
    unsigned int nd = env().getNd();

    if (g->Nd() != nd)
    {
        HADRONS_ERROR(Size, "input field has the wrong number of dimensions");
    }
    fft.FFT_dim_mask(out, in, maskc, backward ? FFT::backward : FFT::forward);
}

template <typename Field>
void TFFT<Field>::execute(void)
{
    std::vector<int> mask = strToVec<int>(par().dimMask);

    if (mask.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "dimension mask has the wrong number of components");
    }
    LOG(Message) << "Dimension mask: " << mask << std::endl;
    LOG(Message) << "     Direction: " << (par().backward ? "backward" : "forward") << std::endl;
    if (isVector_)
    {
        auto &field = envGet(std::vector<Field>, par().field);
        auto &out = envGet(std::vector<Field>, getName());

        LOG(Message) << "Input '" << par().field << "' is a vector of " << field.size() 
                     << " field(s)" << std::endl;
        for (unsigned int i = 0; i < field.size(); ++i)
        {
            LOG(Message) << "Performing FFT on component " << i << std::endl;
            doFft(out[i], field[i], mask, par().backward);
        }
    }
    else
    {
        auto &field = envGet(Field, par().field);
        auto &out = envGet(Field, getName());
        LOG(Message) << "Performing FFT on field '" << par().field << "'" << std::endl;
        doFft(out, field, mask, par().backward);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_FFT_hpp_
