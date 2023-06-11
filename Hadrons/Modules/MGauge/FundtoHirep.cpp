/*
 * FundtoHirep.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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

#include <Hadrons/Modules/MGauge/FundtoHirep.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

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
    envCreateLat(Rep::LatticeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <class Rep>
void TFundtoHirep<Rep>::execute(void)
{
    LOG(Message) << "Transforming Representation" << std::endl;

    auto &U    = envGet(LatticeGaugeField, par().gaugeconf);
    auto &URep = envGet(Rep::LatticeField, getName());

    Rep TargetRepresentation(U._grid);
    TargetRepresentation.update_representation(U);
    URep = TargetRepresentation.U;
}
