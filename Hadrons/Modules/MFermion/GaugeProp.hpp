/*
 * GaugeProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Guido Cossu <guido.cossu@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: pretidav <david.preti@csic.es>
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

#ifndef Hadrons_MFermion_GaugeProp_hpp_
#define Hadrons_MFermion_GaugeProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugeProp                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class GaugePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugePropPar,
                                    std::string, source,
                                    std::string, solver);
};

template <typename FImpl>
class TGaugeProp: public Module<GaugePropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TGaugeProp(const std::string name);
    // destructor
    virtual ~TGaugeProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void solvePropagator(std::vector<PropagatorField *> &prop, 
                         std::vector<PropagatorField *> &propPhysical,
                         const std::vector<PropagatorField *> &fullSrc);
private:
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(GaugeProp, TGaugeProp<FIMPL>, MFermion);
MODULE_REGISTER_TMP(ZGaugeProp, TGaugeProp<ZFIMPL>, MFermion);

/******************************************************************************
 *                      TGaugeProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaugeProp<FImpl>::TGaugeProp(const std::string name)
: Module<GaugePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugeProp<FImpl>::setup(void)
{
    unsigned int sourceSize;

    Ls_ = env().getObjectLs(par().solver);
    // source is a single propagator
    if (envHasType(PropagatorField, par().source))
    {
        envCreateLat(PropagatorField, getName());
        if (Ls_ > 1)
        {
            envCreateLat(PropagatorField, getName() + "_5d", Ls_);
        }
        sourceSize = Ns*FImpl::Dimension;
    }
    // source is a vector of propagators
    else if (envHasType(std::vector<PropagatorField>, par().source))
    {
        auto &src = envGet(std::vector<PropagatorField>, par().source);

        envCreate(std::vector<PropagatorField>, getName(), 1, src.size(),
                  envGetGrid(PropagatorField));
        if (Ls_ > 1)
        {
            envCreate(std::vector<PropagatorField>, getName() + "_5d", Ls_,
                      src.size(), envGetGrid(PropagatorField, Ls_));
        }
        sourceSize = src.size()*Ns*FImpl::Dimension;
    }
    // source is a vector of pointer of propagators
    else if (envHasType(std::vector<PropagatorField *>, par().source))
    {
        auto &src = envGet(std::vector<PropagatorField *>, par().source);

        envCreate(std::vector<PropagatorField>, getName(), 1, src.size(),
                  envGetGrid(PropagatorField));
        if (Ls_ > 1)
        {
            envCreate(std::vector<PropagatorField>, getName() + "_5d", Ls_,
                      src.size(), envGetGrid(PropagatorField, Ls_));
        }
        sourceSize = src.size()*Ns*FImpl::Dimension;
    }
    else
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().source 
                          + "' has an incompatible type ("
                          + env().getObjectType(par().source)
                          + ")", env().getObjectAddress(par().source))
    }
    envTmpLat(FermionField, "tmp");
    if (Ls_ > 1)
    {
        envTmp(std::vector<FermionField>, "source", Ls_, sourceSize,
               envGetGrid(FermionField, Ls_));
        envTmp(std::vector<FermionField>, "sol", Ls_, sourceSize,
               envGetGrid(FermionField, Ls_));
    }
    else
    {
        envTmp(std::vector<FermionField>, "source", 1, sourceSize,
               envGetGrid(FermionField));
        envTmp(std::vector<FermionField>, "sol", 1, sourceSize,
               envGetGrid(FermionField));
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugeProp<FImpl>::solvePropagator(std::vector<PropagatorField *> &prop, 
                                        std::vector<PropagatorField *> &propPhysical,
                                        const std::vector<PropagatorField *> &fullSrc)
{
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();
    unsigned int j = 0;
    
    envGetTmp(std::vector<FermionField>, source);
    envGetTmp(std::vector<FermionField>, sol);
    envGetTmp(FermionField, tmp);

    LOG(Message) << "Import sources" << std::endl;
    startTimer("Import sources");
    for (unsigned int i = 0; i < fullSrc.size(); ++i)
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        // 4D sources
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
               PropToFerm<FImpl>(source[j], *(fullSrc[i]), s, c);
            }
            else
            {
                PropToFerm<FImpl>(tmp, *(fullSrc[i]), s, c);
                mat.ImportPhysicalFermionSource(tmp, source[j]);
            }
        }
        // 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm<FImpl>(source[j], *(fullSrc[i]), s, c);
            }
        }
        j++;
    }
    stopTimer("Import sources");
    LOG(Message) << "Solve" << std::endl;
    startTimer("Solver");
    for (auto &s: sol)
    {
        s = Zero();
    }
    solver(sol, source);
    stopTimer("Solver");
    LOG(Message) << "Export solutions" << std::endl;
    startTimer("Export solutions");
    j = 0;
    for (unsigned int i = 0; i < fullSrc.size(); ++i)
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        FermToProp<FImpl>(*(prop[i]), sol[j], s, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            mat.ExportPhysicalFermionSolution(sol[j], tmp);
            FermToProp<FImpl>(*(propPhysical[i]), tmp, s, c);
        }
        j++;
    }
    stopTimer("Export solutions");
}

template <typename FImpl>
void TGaugeProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    std::string propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    std::vector<PropagatorField *> propPt, physPropPt, srcPt;

    // source is a single propagator
    if (envHasType(PropagatorField, par().source))
    {
        auto &prop         = envGet(PropagatorField, propName);
        auto &propPhysical = envGet(PropagatorField, getName());
        auto &fullSrc      = envGet(PropagatorField, par().source);

        LOG(Message) << "Using source '" << par().source << "'" << std::endl;
        propPt.push_back(&prop);
        physPropPt.push_back(&propPhysical);
        srcPt.push_back(&fullSrc);
    }
    // source is a vector of propagators
    else if (envHasType(std::vector<PropagatorField>, par().source))
    {
        auto &prop         = envGet(std::vector<PropagatorField>, propName);
        auto &propPhysical = envGet(std::vector<PropagatorField>, getName());
        auto &fullSrc      = envGet(std::vector<PropagatorField>, par().source);

        LOG(Message) << "Using source vector '" << par().source << "'" << std::endl;
        for (unsigned int i = 0; i < fullSrc.size(); ++i)
        {
            propPt.push_back(&(prop[i]));
            physPropPt.push_back(&(propPhysical[i]));
            srcPt.push_back(&(fullSrc[i]));
        }
    }
    // source is a vector of pointer of propagators
    else if (envHasType(std::vector<PropagatorField *>, par().source))
    {
        auto &prop         = envGet(std::vector<PropagatorField>, propName);
        auto &propPhysical = envGet(std::vector<PropagatorField>, getName());
        auto &fullSrcPt    = envGet(std::vector<PropagatorField *>, par().source);

        LOG(Message) << "Using source reference vector '" << par().source << "'" << std::endl;
        for (unsigned int i = 0; i < fullSrcPt.size(); ++i)
        {
            propPt.push_back(&(prop[i]));
            physPropPt.push_back(&(propPhysical[i]));
        }
        srcPt = fullSrcPt;
    }
    solvePropagator(propPt, physPropPt, srcPt);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugeProp_hpp_
