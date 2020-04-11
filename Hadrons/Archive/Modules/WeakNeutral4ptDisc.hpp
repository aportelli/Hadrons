/*
 * WeakNeutral4ptDisc.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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

#ifndef Hadrons_MContraction_WeakNeutral4ptDisc_hpp_
#define Hadrons_MContraction_WeakNeutral4ptDisc_hpp_

#include <Hadrons/Modules/MContraction/WeakHamiltonian.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         WeakNeutral4ptDisc                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

enum
{
    neut_disc_1_diag = 0,
    neut_disc_2_diag = 1,
    n_neut_disc_diag = 2
};

// Neutral 4pt disconnected subdiagram contractions.
#define MAKE_DISC_MESON(Q_1, Q_2, gamma) (Q_1*adj(Q_2)*g5*gamma)
#define MAKE_DISC_LOOP(Q_LOOP, gamma) (Q_LOOP*gamma)
#define MAKE_DISC_CURR(Q_c, gamma) (trace(Q_c*gamma))

MAKE_WEAK_MODULE(WeakNeutral4ptDisc)

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakNeutral4ptDisc_hpp_
