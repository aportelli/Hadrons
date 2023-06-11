/*
 * EigenPackLCCompress.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MUtilities_EigenPackLCCompress_hpp_
#define Hadrons_MUtilities_EigenPackLCCompress_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  Local coherence eigenvector compressor                    *
 *****************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class EigenPackLCCompressPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(EigenPackLCCompressPar,
                                    std::string,      epack,
                                    std::string,      blockSize,
                                    unsigned int,     coarseSize,
                                    unsigned int,     Ls,
                                    std::string,      output,
                                    bool,             multiFile);
};

template <typename FImpl, int nBasis, typename FImplIo = FImpl>
class TEigenPackLCCompress: public Module<EigenPackLCCompressPar>
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
    TEigenPackLCCompress(const std::string name);
    // destructor
    virtual ~TEigenPackLCCompress(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(EigenPackLCCompress    , ARG(TEigenPackLCCompress<FIMPL , HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCCompress250 , ARG(TEigenPackLCCompress<FIMPL , 250>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCCompress400 , ARG(TEigenPackLCCompress<FIMPL , 400>), MUtilities);

MODULE_REGISTER_TMP(EigenPackLCCompressF   , ARG(TEigenPackLCCompress<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCCompress250F, ARG(TEigenPackLCCompress<FIMPLF, 250>), MUtilities);
MODULE_REGISTER_TMP(EigenPackLCCompress400F, ARG(TEigenPackLCCompress<FIMPLF, 400>), MUtilities);



/******************************************************************************
 *                 TEigenPackLCCompress implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
TEigenPackLCCompress<FImpl, nBasis, FImplIo>::TEigenPackLCCompress(const std::string name)
: Module<EigenPackLCCompressPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TEigenPackLCCompress<FImpl, nBasis, FImplIo>::getInput(void)
{
    std::vector<std::string> in = {par().epack};
    
    return in;
}

template <typename FImpl, int nBasis, typename FImplIo>
std::vector<std::string> TEigenPackLCCompress<FImpl, nBasis, FImplIo>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TEigenPackLCCompress<FImpl, nBasis, FImplIo>::setup(void)
{
    GridBase *gridIo = nullptr, *gridCoarseIo = nullptr;

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

    envCreate(CoarsePack, getName(), par().Ls,
                     nBasis, par().coarseSize, envGetRbGrid(Field, par().Ls), gridCoarse,
                     gridIo, gridCoarseIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis, typename FImplIo>
void TEigenPackLCCompress<FImpl, nBasis, FImplIo>::execute(void)
{
    auto &coarsePack   = envGet(CoarsePack, getName());
    auto &finePack     = envGet(BasePack, par().epack);

    unsigned int sizeFine = nBasis;
    unsigned int sizeCoarse = par().coarseSize;

    coarsePack.record = finePack.record;

    LOG(Message) << "Copying lowest " << sizeFine << " fine vectors as basis" << std::endl;
    for (unsigned int i=0; i<sizeFine; i++)
    {
        coarsePack.eval[i] = finePack.eval[i];
        coarsePack.evec[i] = finePack.evec[i];
    }

    if (!par().output.empty())
    {
        LOG(Message) << "Write " << sizeFine << " fine basis vectors" << std::endl;
        coarsePack.writeFine(par().output, par().multiFile, vm().getTrajectory());
    }

    auto blockSize = strToVec<int>(par().blockSize);
    GridBase *gridCoarse = envGetCoarseGrid(CoarseField, blockSize, par().Ls);

    LOG(Message) <<"Orthogonalising fine basis"<<std::endl;
    Lattice<typename FImpl::SiteComplex> innerProduct(gridCoarse);
    LOG(Message) <<" Block Gramm-Schmidt pass 1"<<std::endl;
    blockOrthonormalize(innerProduct,coarsePack.evec);
    LOG(Message) <<" Block Gramm-Schmidt pass 2"<<std::endl;
    blockOrthonormalize(innerProduct,coarsePack.evec);

    LOG(Message) << "Projecting " << sizeCoarse << " coarse eigenvectors" << std::endl;
    for (unsigned int i=0; i<finePack.evec.size(); i++)
    {
        LOG(Message) << "evec " << i << std::endl;
        blockProject(coarsePack.evecCoarse[i], finePack.evec[i], coarsePack.evec);
        coarsePack.evalCoarse[i] = finePack.eval[i];
    }

    if (!par().output.empty())
    {
        LOG(Message) << "Write " << sizeCoarse << " coarse vectors" << std::endl;
        coarsePack.writeCoarse(par().output, par().multiFile, vm().getTrajectory());
    }


}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_EigenPackLCCompress_hpp_
