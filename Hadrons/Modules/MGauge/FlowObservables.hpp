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
                                    unsigned int, Nstep);
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

    Usmear = U;
    std::vector<std::pair<RealD,RealD>> plaquette_density;
    std::vector<std::pair<RealD,RealD>> clover_density;
    std::vector<std::pair<RealD,RealD>> topq;

    plaquette_density.resize(0);
    clover_density.resize(0);
    topq.resize(0);

    WilsonFlow<PeriodicGimplR> wflow(par().epsilon, par().Nstep, 1);
    wflow.resetActions();

    wflow.addMeasurement(1, [&wflow, &clover_density, &plaquette_density](int step, RealD t, const LatticeGaugeField &U){
        RealD tmp_plaquette = wflow.energyDensityPlaquette(t,U);
        LOG(Debug) << "Plaquette energy density for step " << step << ", t=" << t << ", E=" << tmp_plaquette << std::endl;
        plaquette_density.push_back( {t, tmp_plaquette} );

        RealD tmp_clover = wflow.energyDensityCloverleaf(t,U);
        LOG(Debug) << "Cloverleaf energy density for step " << step << ", t=" << t << ", E=" << tmp_clover << std::endl;
        clover_density.push_back( {t, tmp_clover} );
    });
    wflow.addMeasurement(1, [&wflow, &topq](int step, RealD t, const LatticeGaugeField &U){
        RealD tmp_topq = WilsonLoops<GImpl>::TopologicalCharge(U);
        LOG(Debug) << "Topological charge for step " << step << ", t=" << t << ", Q=" << tmp_topq << std::endl;
        topq.push_back({t, tmp_topq});
    });

    wflow.smear(Usmear, U);

    Real qtop = WilsonLoops<GImpl>::TopologicalCharge(Usmear);

    LOG(Message) << "Old method Topological charge: " << qtop << std::endl;
    LOG(Message) << "Topological charge " << topq.back().second << " at flow time " << topq.back().first << std::endl;

    LOG(Message) << "Energy density:" << std::endl;
    for (auto en: clover_density)
        LOG(Message) << en << std::endl;

    LOG(Message) << "Topological charge:" << std::endl;
    for (auto en: topq)
        LOG(Message) << en << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_FlowObservables_hpp_
