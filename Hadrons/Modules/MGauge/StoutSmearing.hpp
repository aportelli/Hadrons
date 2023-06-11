/*
 * StoutSmearing.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
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
#ifndef Hadrons_MGauge_StoutSmearing_hpp_
#define Hadrons_MGauge_StoutSmearing_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            Stout smearing                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class StoutSmearingPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StoutSmearingPar,
                                    std::string, gauge,
                                    unsigned int, steps,
                                    std::string, orthogDim,
                                    double, rho);
};

template <typename GImpl>
class TStoutSmearing: public Module<StoutSmearingPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TStoutSmearing(const std::string name);
    // destructor
    virtual ~TStoutSmearing(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StoutSmearing, TStoutSmearing<GIMPL>, MGauge);

/******************************************************************************
 *                     TStoutSmearing implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TStoutSmearing<GImpl>::TStoutSmearing(const std::string name)
: Module<StoutSmearingPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TStoutSmearing<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl>
std::vector<std::string> TStoutSmearing<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TStoutSmearing<GImpl>::setup(void)
{
    envCreateLat(GaugeField, getName());
    envTmpLat(GaugeField, "buf");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TStoutSmearing<GImpl>::execute(void)
{
    LOG(Message) << "Smearing '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps > 1) ? "s" : "") 
                 << " of stout smearing and rho= " << par().rho << std::endl;

    int oDim = -1;
    if(!par().orthogDim.empty())
    {
        LOG(Message) << "Only smearing orthogonally to dimension " << par().orthogDim << std::endl;
        oDim=std::stoi(par().orthogDim);
    }
    Smear_Stout<GImpl> smearer(par().rho, oDim);
    auto               &U    = envGet(GaugeField, par().gauge);
    auto               &Usmr = envGet(GaugeField, getName());

    envGetTmp(GaugeField, buf);
    buf = U;
    LOG(Message) << "plaquette= " << WilsonLoops<GImpl>::avgPlaquette(U)
                 << std::endl;
    // This is only here to return the original field in case steps=0. Not sure whether it's needed.
    Usmr = U;
    for (unsigned int n = 0; n < par().steps; ++n)
    {
        smearer.smear(Usmr, buf);
        buf = Usmr;
        LOG(Message) << "plaquette= " << WilsonLoops<GImpl>::avgPlaquette(Usmr)
                     << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_StoutSmearing_hpp_
