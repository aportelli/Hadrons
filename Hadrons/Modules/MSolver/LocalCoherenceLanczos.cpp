/*
 * LocalCoherenceLanczos.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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
#include <Hadrons/Modules/MSolver/LocalCoherenceLanczos.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<FIMPL,600>;
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>;
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<FIMPL,600, FIMPLF>;
template class Grid::Hadrons::MSolver::TLocalCoherenceLanczos<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, ZFIMPLF>;
#endif
