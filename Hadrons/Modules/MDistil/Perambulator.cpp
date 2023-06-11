/*
 * Perambulator.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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

#include <Hadrons/Modules/MDistil/Perambulator.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MDistil;

template class HADRONS_NAMESPACE::MDistil::TPerambulator<FIMPL>;
template class HADRONS_NAMESPACE::MDistil::TPerambulator<ZFIMPL>;

BEGIN_HADRONS_NAMESPACE

// Global constants for distillation

BEGIN_MODULE_NAMESPACE(MDistil)

const std::string                PerambTensor::Name__{"Perambulator"};
const std::array<std::string, 6> PerambTensor::DefaultIndexNames__{"nT", "nVec", "LI", "nNoise", "nT_inv", "SI"};

const std::string                PerambIndexTensor::Name__{"Perambulator"};
const std::array<std::string, 5> PerambIndexTensor::DefaultIndexNames__{"nT", "nVec", "nDL", "nNoise", "nDS"};

const std::string                TimesliceEvals::Name__{"TimesliceEigenValues"};
const std::array<std::string, 2> TimesliceEvals::DefaultIndexNames__{"nT", "nVec"};

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
