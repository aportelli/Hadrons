/*
 * LoadCoarseEigenPack.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_LoadCoarseEigenPack_hpp_
#define Hadrons_MIO_LoadCoarseEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Load local coherence eigen vectors/values package             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCoarseEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCoarseEigenPackPar,
                                    std::string,  filestem,
                                    bool,         multiFile,
                                    unsigned int, sizeFine,
                                    unsigned int, sizeCoarse,
                                    bool,         redBlack,
                                    unsigned int, Ls,
                                    std::string,  blockSize,
                                    bool,         orthogonalise,
                                    std::string,  gaugeXform);
};

template <typename Pack, typename GImpl>
class TLoadCoarseEigenPack: public Module<LoadCoarseEigenPackPar>
{
public:
    typedef typename Pack::Field                Field;
    typedef typename Pack::FieldIo              FieldIo;
    typedef typename Pack::CoarseField          CoarseField;
    typedef typename Pack::CoarseFieldIo        CoarseFieldIo;
    typedef CoarseEigenPack<Field, CoarseField, FieldIo, CoarseFieldIo> BasePack;
    template <typename vtype> 
    using iImplScalar = iScalar<iScalar<iScalar<vtype>>>;
    typedef iImplScalar<typename Pack::Field::vector_type> SiteComplex;

    GAUGE_TYPE_ALIASES(GImpl, );
    typedef typename GImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TLoadCoarseEigenPack(const std::string name);
    // destructor
    virtual ~TLoadCoarseEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack   , ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>, GIMPL>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack250, ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, 250>, GIMPL>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack400, ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL, 400>, GIMPL>), MIO);

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPackF   , ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>, GIMPLF>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack250F, ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF, 250>, GIMPLF>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPack400F, ARG(TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF, 400>, GIMPLF>), MIO);


/******************************************************************************
 *                 TLoadCoarseEigenPack implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
TLoadCoarseEigenPack<Pack,GImpl>::TLoadCoarseEigenPack(const std::string name)
: Module<LoadCoarseEigenPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
std::vector<std::string> TLoadCoarseEigenPack<Pack,GImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().gaugeXform.empty())
    {
        in = {par().gaugeXform};
    }
    
    return in;
}

template <typename Pack, typename GImpl>
std::vector<std::string> TLoadCoarseEigenPack<Pack,GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadCoarseEigenPack<Pack,GImpl>::setup(void)
{
    GridBase     *gridIo = nullptr, *gridCoarseIo = nullptr;

    auto blockSize = strToVec<int>(par().blockSize);

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = envGetRbGrid(FieldIo, par().Ls);
    }
    if (typeHash<CoarseField>() != typeHash<CoarseFieldIo>())
    {
        gridCoarseIo = envGetCoarseGrid(CoarseFieldIo, blockSize, par().Ls);
    }
    envCreateDerived(BasePack, Pack, getName(), par().Ls, par().sizeFine,
                     par().sizeCoarse, envGetRbGrid(Field, par().Ls), 
                     envGetCoarseGrid(CoarseField, blockSize, par().Ls),
                     gridIo, gridCoarseIo);

    GridBase* grid   = envGetGrid(Field, par().Ls);
    GridBase* gridRb = envGetRbGrid(Field, par().Ls);
    if (!par().gaugeXform.empty())
    {
        envTmp(GaugeMat,    "tmpXform", par().Ls, grid);
        if (par().redBlack)
        {
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, gridRb);
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadCoarseEigenPack<Pack,GImpl>::execute(void)
{
    auto blockSize = strToVec<int>(par().blockSize);

    auto                 cg     = envGetCoarseGrid(CoarseField, blockSize, par().Ls);
    auto                 &epack = envGetDerived(BasePack, Pack, getName());
    Lattice<SiteComplex> dummy(cg);

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());

    if (!par().gaugeXform.empty())
    {
        LOG(Message) << "Applying gauge transformation to fine eigenvectors " << getName()
                     << " using " << par().gaugeXform << std::endl;

        auto &xform = envGet(GaugeMat, par().gaugeXform);
        epack.gaugeTransform(xform);
    }

    if (par().orthogonalise) {
        LOG(Message) << "Block Gramm-Schmidt pass 1"<< std::endl;
        blockOrthogonalise(dummy, epack.evec);
        LOG(Message) << "Block Gramm-Schmidt pass 2"<< std::endl;
        blockOrthogonalise(dummy, epack.evec);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCoarseEigenPack_hpp_
