/*
 * FundtoHirep.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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

#ifndef Hadrons_MGauge_FundtoHirep_hpp_
#define Hadrons_MGauge_FundtoHirep_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         Load a NERSC configuration                         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class FundtoHirepPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FundtoHirepPar,
                                    std::string, gaugeconf);
};

template <class Rep>
class TFundtoHirep: public Module<FundtoHirepPar>
{
public:
    // constructor
    TFundtoHirep(const std::string name);
    // destructor
    virtual ~TFundtoHirep(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    void setup(void);
    // execution
    void execute(void);
};

MODULE_REGISTER_TMP(FundtoAdjoint,   TFundtoHirep<AdjointRepresentation>, MGauge);
MODULE_REGISTER_TMP(FundtoTwoIndexSym, TFundtoHirep<TwoIndexSymmetricRepresentation>, MGauge);
MODULE_REGISTER_TMP(FundtoTwoIndexAsym, TFundtoHirep<TwoIndexAntiSymmetricRepresentation>, MGauge);

// constructor /////////////////////////////////////////////////////////////////
template <class Rep>
TFundtoHirep<Rep>::TFundtoHirep(const std::string name)
: Module<FundtoHirepPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <class Rep>
std::vector<std::string> TFundtoHirep<Rep>::getInput(void)
{
    std::vector<std::string> in = {par().gaugeconf};

    return in;
}

template <class Rep>
std::vector<std::string> TFundtoHirep<Rep>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Rep>
void TFundtoHirep<Rep>::setup(void)
{
    envCreateLatNS(Rep::LatticeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <class Rep>
void TFundtoHirep<Rep>::execute(void)
{
    LOG(Message) << "Transforming Representation" << std::endl;

    auto &U    = envGet(LatticeGaugeField, par().gaugeconf);
    auto &URep = envGet(typename Rep::LatticeField, getName());

    Rep TargetRepresentation(U.Grid());
    TargetRepresentation.update_representation(U);
    URep = TargetRepresentation.U;
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_FundtoHirep_hpp_
