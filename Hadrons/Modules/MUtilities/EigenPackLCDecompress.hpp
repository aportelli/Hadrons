/*
 * EigenPackLCDecompress.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MUtilities_EigenPackLCDecompress_hpp_
#define Hadrons_MUtilities_EigenPackLCDecompress_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Local coherence eigenvector decompressor                   *
 *****************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EigenPackLCDecompressPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackLCDecompressPar,
                                    std::string,      epack,
                                    std::string,      blockSize,
                                    unsigned int,     coarseSize,
                                    unsigned int,     Ls,
                                    std::string,      output,
                                    bool,             multiFile);
};

template <typename FImpl, int nBasis, typename FImplIo = FImpl>
class TEigenPackLCDecompress: public Module<EigenPackLCDecompressPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef BaseFermionEigenPack<FImpl>                    BasePack;
    typedef CoarseFermionEigenPack<FImpl, nBasis, FImplIo> CoarsePack;
    typedef typename CoarsePack::Field                     Field;
    typedef typename CoarsePack::FieldIo                   FieldIo;
    typedef typename CoarsePack::CoarseField               CoarseField;
    typedef typename CoarsePack::CoarseFieldIo             CoarseFieldIo;

    typedef FermionEigenPack<FImpl>                        FinePack;
public:
    // constructor
    TEigenPackLCDecompress(const std::string name);
    // destructor
    virtual ~TEigenPackLCDecompress(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(EigenPackLCDecompress    , ARG(TEigenPackLCDecompress<FIMPL , HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCDecompress250 , ARG(TEigenPackLCDecompress<FIMPL , 250>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCDecompress400 , ARG(TEigenPackLCDecompress<FIMPL , 400>), MUtilities);

MODULE_REGISTER_TMP(EigenPackLCDecompressF   , ARG(TEigenPackLCDecompress<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCDecompress250F, ARG(TEigenPackLCDecompress<FIMPLF, 250>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCDecompress400F, ARG(TEigenPackLCDecompress<FIMPLF, 400>), MUtilities);


/******************************************************************************
 *                 TEigenPackLCDecompress implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
TEigenPackLCDecompress<FImpl, nBasis, FImplIo>::TEigenPackLCDecompress(const std::string name)
: Module<EigenPackLCDecompressPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TEigenPackLCDecompress<FImpl, nBasis, FImplIo>::getInput(void)
{
    std::vector<std::string> in = {par().epack};
    
    return in;
}

template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TEigenPackLCDecompress<FImpl, nBasis, FImplIo>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TEigenPackLCDecompress<FImpl, nBasis, FImplIo>::setup(void)
{    
    GridBase *gridIo = nullptr;
    [[maybe_unused]] GridBase *gridCoarseIo = nullptr;

    auto blockSize = strToVec<int>(par().blockSize);

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = envGetRbGrid(FieldIo, par().Ls);
    }
    if (typeHash<CoarseField>() != typeHash<CoarseFieldIo>())
    {
        gridCoarseIo = envGetCoarseGrid(CoarseFieldIo, blockSize, par().Ls);
    }

    GridBase* gridCoarse  = envGetCoarseGrid(CoarseField, blockSize, par().Ls);
    GridBase* grid        = getGrid<Field>(true, par().Ls);

    envCreateDerived(BasePack, FinePack, getName(), par().Ls, par().coarseSize, grid, gridIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TEigenPackLCDecompress<FImpl, nBasis, FImplIo>::execute(void)
{
    auto &coarsePack = envGet(CoarsePack, par().epack);
    auto &finePack   = envGetDerived(BasePack, FinePack, getName());

    finePack.record = coarsePack.record;

    auto blockSize = strToVec<int>(par().blockSize);
    GridBase *gridCoarse = envGetCoarseGrid(CoarseField, blockSize, par().Ls);

    LOG(Message) << "Promoting from coarse basis" << std::endl;
    for (unsigned int i=0; i<finePack.evec.size(); i++)
    {
        LOG(Message) << "evec " << i << std::endl;
        blockPromote(coarsePack.evecCoarse[i], finePack.evec[i], coarsePack.evec);
        finePack.eval[i] = coarsePack.evalCoarse[i];
    }

    if (!par().output.empty())
    {
        LOG(Message) << "Write " << par().coarseSize << " decompressed vectors" << std::endl;
        finePack.write(par().output, par().multiFile, vm().getTrajectory());
    }


}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_EigenPackLCDecompress_hpp_
