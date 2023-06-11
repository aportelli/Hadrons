/*
 * LeptonImpl.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#include <Hadrons/Global.hpp>
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/implementation/WilsonFermion5DImplementation.h>
#include <Grid/qcd/action/fermion/implementation/CayleyFermion5DImplementation.h>
#include <Grid/qcd/action/fermion/implementation/CayleyFermion5Dcache.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsHandImplementation.h>

#ifndef AVX512
#ifndef QPX
#ifndef A64FX
#ifndef A64FXFIXEDSIZE
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsAsmImplementation.h>
#endif
#endif
#endif
#endif

NAMESPACE_BEGIN(Grid);

template class CayleyFermion5D<LeptonWilsonImplD>; 
template class CayleyFermion5D<LeptonWilsonImplF>; 
template class WilsonFermion5D<LeptonWilsonImplD>; 
template class WilsonFermion5D<LeptonWilsonImplF>; 
template class WilsonKernels<LeptonWilsonImplD>;
template class WilsonKernels<LeptonWilsonImplF>;

NAMESPACE_END(Grid);
