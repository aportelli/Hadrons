/*
 * LatticeUtilities.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_LatticeUtilities_hpp_
#define Hadrons_LatticeUtilities_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

// fill a (N+1)-dimension lattice by copying slices of a N-dimension one
// possibly not super optimal
template <typename Field>
void sliceFill(Field &out, const Field &in, const unsigned int dir = 0,
               const int cbComponent = Even)
{
    GridBase *hg = out.Grid(), *g = in.Grid(), *tmpGrid;
    std::shared_ptr<GridBase> tmpGridPt;

    // same grid? simple just copy
    if (hg == g)
    {
        out = in;
    }
    // if not several cases are possible
    else
    {
        unsigned int extraDim;
        bool         cb = false, hcb = false;

        // check if the in field is checkerboarded (not supported)
        for (unsigned int mu = 0; mu < g->Dimensions(); ++mu)
        {
            cb = cb or (g->CheckerBoarded(mu) != 0);
        }
        for (unsigned int mu = 0; mu < hg->Dimensions(); ++mu)
        {
            hcb = hcb or (hg->CheckerBoarded(mu) != 0);
        }
        if (cb)
        {
            HADRONS_ERROR(Implementation, "sliceFill not implemented for checkerboarded input");
        }

        // check is the ouput dimensions are 1+
        if (hg->Dimensions() == g->Dimensions() + 1)
        {
            auto hdim = hg->FullDimensions();

            extraDim = hdim[dir];
            tmpGridPt.reset(new GridCartesian(hdim, hg->_simd_layout, hg->_processors));
            tmpGrid = tmpGridPt.get();
        }
        else
        {
            extraDim = 1;
            tmpGrid  = g;
        }

        // this temporary has the same dimensions than the output
        Field tmp(tmpGrid);

        // if output dimensions are 1+ copy across the extra dimension
        if (extraDim != 1)
        {
            for (unsigned int s = 0; s < extraDim; ++s)
            {
                InsertSlice(in, tmp, s, dir);
            }
        }
        // else just copy
        else
        {
            tmp = in;
        }

        // finally, if checkerboarded pick the right component
        if (hcb)
        {
            pickCheckerboard(cbComponent, out, tmp);
        }
        // else just copy
        else
        {
            out = tmp;
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_LatticeUtilities_hpp_
