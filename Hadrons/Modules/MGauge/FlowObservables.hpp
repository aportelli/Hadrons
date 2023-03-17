/*
 * FlowObservables.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
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

#ifndef Hadrons_MGauge_FlowObservables_hpp_
#define Hadrons_MGauge_FlowObservables_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/smearing/WilsonFlow.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         FlowObservables                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class FlowObservablesPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FlowObservablesPar,
                                    std::string, gauge,
                                    Real, epsilon,
                                    Real, maxTau);
};

template <typename GImpl>
class TFlowObservables: public Module<FlowObservablesPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TFlowObservables(const std::string name);
    // destructor
    virtual ~TFlowObservables(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FlowObservables, TFlowObservables<GIMPL>, MGauge);

/******************************************************************************
 *                 TFlowObservables implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TFlowObservables<GImpl>::TFlowObservables(const std::string name)
: Module<FlowObservablesPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TFlowObservables<GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    return in;
}

template <typename GImpl>
std::vector<std::string> TFlowObservables<GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TFlowObservables<GImpl>::setup(void)
{
    envTmpLat(GaugeField, "Usmear");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TFlowObservables<GImpl>::execute(void)
{

    auto &U = envGet(GaugeField, par().gauge);

    envGetTmp(GaugeField, Usmear);

    WilsonFlow<PeriodicGimplR> WF(0, par().epsilon, 1);
    WF.smear_adaptive(Usmear, U, par().maxTau);

    Real qtop = WilsonLoops<GImpl>::TopologicalCharge(Usmear);

    LOG(Message) << "Topological charge: " << qtop << std::endl;

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_FlowObservables_hpp_
