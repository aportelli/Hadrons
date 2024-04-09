#ifndef Hadrons_MScalar_Product_hpp_
#define Hadrons_MScalar_Product_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                             Field product                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class ProductPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                    std::vector<std::string>, fields);
};

template <typename SImpl>
class TProduct: public Module<ProductPar>
{
public:
    BASIC_TYPE_ALIASES(SImpl,);
public:
    // constructor
    TProduct(const std::string name);
    // destructor
    virtual ~TProduct(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Product, TProduct<SIMPL>, MScalar);

/******************************************************************************
 *                         TProduct implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TProduct<SImpl>::TProduct(const std::string name)
: Module<ProductPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TProduct<SImpl>::getInput(void)
{
    std::vector<std::string> in = par().fields;
    
    return in;
}

template <typename SImpl>
std::vector<std::string> TProduct<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TProduct<SImpl>::setup(void)
{
    envCreateLat(ScalarField, getName());
    if (par().fields.size() == 0)
    {
        HADRONS_ERROR(Size, "field vector input empty");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TProduct<SImpl>::execute(void)
{
    auto &result = envGet(ScalarField, getName());
    auto &phi1 = envGet(ScalarField, par().fields[0]);

    result = phi1;
    for (unsigned int j = 0; j < par().fields.size(); ++j)
    {
        auto &phij = envGet(ScalarField, par().fields[j]);
        result *= phij;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_Product_hpp_
