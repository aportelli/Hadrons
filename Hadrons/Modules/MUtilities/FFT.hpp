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
};

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
    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TFFT<Field>::execute(void)
{
    auto &field = envGet(Field, par().field);
    auto &out = envGet(Field, getName());
    GridBase *g = field.Grid();
    unsigned int nd = env().getNd(), nt = env().getDim(Tp);
    FFT fft(dynamic_cast<GridCartesian *>(g));
    std::vector<int> maskv = strToVec<int>(par().dimMask);
    Coordinate mask(maskv);

    if (g->Nd() != nd)
    {
        HADRONS_ERROR(Size, "input field has the wrong number of dimensions");
    }
    if (maskv.size() != nd)
    {
        HADRONS_ERROR(Size, "dimension mask has the wrong number of components");
    }
    LOG(Message) << "Performing FFT on field '" << par().field << "'" << std::endl;
    LOG(Message) << "Dimension mask: " << maskv << std::endl;
    LOG(Message) << "     Direction: " << (par().backward ? "backward" : "forward") << std::endl;
    fft.FFT_dim_mask(out, field, mask, par().backward ? FFT::backward : FFT::forward);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_FFT_hpp_
