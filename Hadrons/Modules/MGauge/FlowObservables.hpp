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
 *                Computation of Wilson flow observables                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class FlowObservablesPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FlowObservablesPar,
                                    std::string, gauge,
                                    Real, epsilon,
                                    unsigned int, Nstep,
                                    std::string, output);
};

template <typename GImpl>
class TFlowObservables: public Module<FlowObservablesPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
    class FlowResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FlowResult,
                                        std::string, description,
                                        std::vector<RealD>, data);
    };
public:
    // constructor
    TFlowObservables(const std::string name);
    // destructor
    virtual ~TFlowObservables(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
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

template <typename GImpl>
std::vector<std::string> TFlowObservables<GImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
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

    std::vector<FlowResult> flow_result;
    flow_result.resize(4);

    flow_result[0].description = "Flow time";
    flow_result[1].description = "Plaquette energy density";
    flow_result[2].description = "Clover energy density";
    flow_result[3].description = "Topological charge";

    WilsonFlow<PeriodicGimplR> wflow(par().epsilon, par().Nstep, 1);
    wflow.resetActions();

    wflow.addMeasurement(1, [&wflow, &flow_result](int step, RealD t, const LatticeGaugeField &U){
        flow_result[0].data.push_back(t);
        flow_result[1].data.push_back(wflow.energyDensityPlaquette(t,U));
        flow_result[2].data.push_back(wflow.energyDensityCloverleaf(t,U));
        flow_result[3].data.push_back(WilsonLoops<GImpl>::TopologicalCharge(U));

        LOG(Debug) << "Step " << step << std::endl;
        for (unsigned int i = 0; i < 4; ++i)
        {
            LOG(Debug) << std::setw(25) << flow_result[i].description << " " << flow_result[i].data.back() << std::endl;
        }
        LOG(Debug) << std::endl;
    });

    wflow.smear(Usmear, U);

    for (unsigned int i = 1; i < 4; ++i)
    {
        LOG(Message) << std::setw(25) << flow_result[i].description << " : " << flow_result[i].data.back() << std::endl;
    }

    saveResult(par().output, "FlowObservables", flow_result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_FlowObservables_hpp_
