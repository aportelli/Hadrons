#ifndef Hadrons_MCovariantFourier_Decompress_hpp_
#define Hadrons_MCovariantFourier_Decompress_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MCovariantFourier/Spectrum.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Decompress                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MCovariantFourier)

class DecompressPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DecompressPar,
                                    std::string, basis,
                                    std::string, spectrum,
                                    unsigned int, Ls,
                                    unsigned int, size);
};

template <typename FImpl>
class TDecompress: public Module<DecompressPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TDecompress(const std::string name);
    // destructor
    virtual ~TDecompress(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Decompress, TDecompress<FIMPL>, MCovariantFourier);

/******************************************************************************
 *                 TDecompress implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDecompress<FImpl>::TDecompress(const std::string name)
: Module<DecompressPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDecompress<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().basis};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDecompress<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDecompress<FImpl>::setup(void)
{
    GridBase *grid = getGrid<FermionField>(par().Ls);

    envCreateDerived(BaseEigenPack<FermionField>, EigenPack<FermionField>, getName(), par().Ls, par().size, grid);
    envTmpLat(ColourVectorField, "tmp", par().Ls);
    if (par().Ls > 1)
    {
        envTmpLat(ColourVectorField, "vec4");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDecompress<FImpl>::execute(void)
{
    ResultReader   reader(resultFilename(par().spectrum));
    SpectrumResult spec;
    auto           &pack  = envGet(BaseEigenPack<FermionField>, getName());
    auto           &basis = envGet(BaseEigenPack<ColourVectorField>, par().basis);

    envGetTmp(ColourVectorField, tmp);
    read(reader, "spectrum", spec);
    pack.eval = spec.vectorEval;
    for (unsigned int i = 0; i < pack.evec.size(); ++i)
    {
        LOG(Message) << "vector " << i << " decompression" << std::endl;
        if (par().Ls == 1)
        {
            for (unsigned int s = 0; s < Ns; ++s)
            {
                LOG(Message) << "spin= " << s << std::endl;
                tmp = Zero();
                for (unsigned int j = 0; j < basis.evec.size(); ++j)
                {
                    tmp += spec.spectrum(i, j, s, 0)*basis.evec[j];
                }
                pokeSpin(pack.evec[i], tmp, s);
            }
        }
        else
        {
            envGetTmp(ColourVectorField, vec4);
            for (unsigned int s = 0; s < Ns; ++s)
            {
                for (unsigned int t = 0; t < par().Ls; ++t)
                {
                    LOG(Message) << "spin= " << s << " / s= " << t << std::endl;
                    vec4 = Zero();

                    for (unsigned int j = 0; j < basis.evec.size(); ++j)
                    {
                        vec4 += spec.spectrum(i, j, s, t)*basis.evec[j];
                    }
                    InsertSlice(vec4, tmp, t, 0);
                }
                pokeSpin(pack.evec[i], tmp, s);
            }
        }
    }
    for (unsigned int i = 0; i < par().size; ++i)
    {
        LOG(Message) << " ";
        for (unsigned int j = 0; j < par().size; ++j)
        {
            std::cout << innerProduct(pack.evec[i], pack.evec[j]) << " ";
        }
        std::cout << std::endl;
    }
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MCovariantFourier_Decompress_hpp_
