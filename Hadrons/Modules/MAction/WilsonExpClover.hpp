/*
 * WilsonExpClover.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Guido Cossu <guido.cossu@ed.ac.uk>
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

#ifndef Hadrons_MAction_WilsonExpClover_hpp_
#define Hadrons_MAction_WilsonExpClover_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Wilson exponential clover quark action                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class WilsonExpCloverPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WilsonExpCloverPar,
                                    std::string, gauge,
                                    double     , mass,
				                    double     , csw_r,
				                    double     , csw_t,
                                    double     , cF,
				                    WilsonAnisotropyCoefficients ,clover_anisotropy,
                                    std::string, boundary,
                                    std::string, twist
				    );
};

template <typename FImpl>
class TWilsonExpClover: public Module<WilsonExpCloverPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWilsonExpClover(const std::string name);
    // destructor
    virtual ~TWilsonExpClover(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonExpClover, TWilsonExpClover<FIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(WilsonExpCloverF, TWilsonExpClover<FIMPLF>, MAction);
#endif

/******************************************************************************
 *                    TWilsonExpClover template implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWilsonExpClover<FImpl>::TWilsonExpClover(const std::string name)
: Module<WilsonExpCloverPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWilsonExpClover<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};

    return in;
}

template <typename FImpl>
std::vector<std::string> TWilsonExpClover<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilsonExpClover<FImpl>::setup(void)
{
    LOG(Message) << "Setting up Wilson exponential clover fermion matrix with m = " << par().mass
                 << " using gauge field '" << par().gauge << "'" << std::endl;
    LOG(Message) << "Clover term csw_r: " << par().csw_r
                 << " csw_t: " << par().csw_t << std::endl;
    LOG(Message) << "Boundary improvement coefficient cF = " << par().cF
                 << std::endl;
                 
    auto &U      = envGet(GaugeField, par().gauge);
    auto &grid   = *envGetGrid(FermionField);
    auto &gridRb = *envGetRbGrid(FermionField);
    typename CompactWilsonCloverFermion<FImpl, CompactExpCloverHelpers<FImpl>>::ImplParams implParams;
    if (!par().boundary.empty())
    {
        implParams.boundary_phases = strToVec<Complex>(par().boundary);
    }
    if (!par().twist.empty())
    {
        implParams.twist_n_2pi_L   = strToVec<Real>(par().twist);
    }
    LOG(Message) << "Fermion boundary conditions: " << implParams.boundary_phases
                 << std::endl;
    LOG(Message) << "Twists: " << implParams.twist_n_2pi_L
                 << std::endl;
    if (implParams.boundary_phases.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of boundary phase");
    }
    if (implParams.twist_n_2pi_L.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "Wrong number of twist");
    }
    envCreateDerived(FMat, CompactWilsonExpClover<FImpl>, getName(), 1, U, grid,
                     gridRb, par().mass, par().csw_r, par().csw_t, par().cF, 
                     par().clover_anisotropy, implParams); 
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWilsonExpClover<FImpl>::execute()
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WilsonExpClover_hpp_
