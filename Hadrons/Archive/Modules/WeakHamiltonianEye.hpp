/*
 * WeakHamiltonianEye.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
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

#ifndef Hadrons_MContraction_WeakHamiltonianEye_hpp_
#define Hadrons_MContraction_WeakHamiltonianEye_hpp_

#include <Hadrons/Modules/MContraction/WeakHamiltonian.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WeakHamiltonianEye                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

enum
{
    S_diag = 0,
    E_diag = 1,
    n_eye_diag = 2
};

// Saucer and Eye subdiagram contractions.
#define MAKE_SE_BODY(Q_1, Q_2, Q_3, gamma) (Q_3*g5*Q_1*adj(Q_2)*g5*gamma)
#define MAKE_SE_LOOP(Q_loop, gamma) (Q_loop*gamma)

MAKE_WEAK_MODULE(WeakHamiltonianEye)

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakHamiltonianEye_hpp_
