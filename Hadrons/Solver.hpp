/*
 * Solver.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_Solver_hpp_
#define Hadrons_Solver_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

template <typename FImpl>
class Solver
{
public:
    typedef typename FImpl::FermionField                      FermionField;
    typedef FermionOperator<FImpl>                            FMat; 
    typedef std::function<void(FermionField &, 
                               const FermionField &)>         SolverFn;
    typedef std::function<void(std::vector<FermionField> &, 
                               const std::vector<FermionField> &)> SolverVFn;
public:
    Solver(SolverFn fn, FMat &mat): mat_(mat), fn_(fn), vfn_(nullptr) {}
    Solver(SolverFn fn, SolverVFn vfn, FMat &mat): mat_(mat), fn_(fn), vfn_(vfn) {}

    void operator()(FermionField &sol, const FermionField &src)
    {
        fn_(sol, src);
    }

    void operator()(std::vector<FermionField> &sol, const std::vector<FermionField> &src)
    {
        if (sol.size() != src.size())
        {
            HADRONS_ERROR(Size, "source and solution vectors size mismatch");
        }
        if (vfn_ == nullptr)
        {
            for (unsigned int i = 0; i < sol.size(); ++i)
            {
                fn_(sol[i], src[i]);
            }
        }
        else
        {
            vfn_(sol, src);
        }
    }

    FMat & getFMat(void)
    {
        return mat_;
    }
private:
    FMat      &mat_;
    SolverFn  fn_;
    SolverVFn vfn_;
};

END_HADRONS_NAMESPACE

#endif // Hadrons_Solver_hpp_
