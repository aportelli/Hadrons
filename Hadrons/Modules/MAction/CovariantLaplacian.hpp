/*
 * CovariantLaplacian.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MAction_CovariantLaplacian_hpp_
#define Hadrons_MAction_CovariantLaplacian_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Covariant Laplacian operator                          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class CovariantLaplacianPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CovariantLaplacianPar,
                                    std::string, gauge,
                                    double,      m2);
};

template <typename Field, typename GImpl>
class CovariantLaplacian: public LinearOperatorBase<Field>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    CovariantLaplacian(const GaugeField &U, const double m2)
    : grid_(U.Grid()), tmp_(U.Grid()), m2_(m2)
    {
        for (unsigned int mu = 0; mu < grid_->Nd(); ++mu)
        {
            U_.push_back(PeekIndex<LorentzIndex>(U, mu)); 
        }
    }

    virtual void OpDiag(const Field &in, Field &out)
    {
        HADRONS_ERROR(Implementation, "CovariantLaplacian method not implemented");
    }

    virtual void OpDir(const Field &in, Field &out, int dir, int disp)
    {
        HADRONS_ERROR(Implementation, "CovariantLaplacian method not implemented");
    }

    virtual void OpDirAll(const Field &in, std::vector<Field> &out)
    {
        HADRONS_ERROR(Implementation, "CovariantLaplacian method not implemented");
    }

    virtual void Op(const Field &in, Field &out)
    {
        unsigned int nd = grid_->Nd();

        out = (m2_ + 2.*nd)*in;
        for (unsigned int mu = 0; mu < nd; ++mu)
        {
            out  -= U_[mu]*Cshift(in, mu, 1);
            tmp_  = adj(U_[mu])*in;
            out  -= Cshift(tmp_, mu, -1);
        }
    }

    virtual void AdjOp(const Field &in, Field &out)
    {
        Op(in, out);
    }

    virtual void HermOpAndNorm(const Field &in, Field &out, RealD &n1, RealD &n2)
    {
        HermOp(in,out);
        Complex dot = innerProduct(in,out); 
        n1 = real(dot);
        n2 = norm2(out);
    }

    virtual void HermOp(const Field &in, Field &out)
    {
        Op(in, out);
    }
private:
    std::vector<GaugeLinkField> U_;
    Field                       tmp_;
    GridBase                    *grid_;
    double                      m2_;
};

template <typename Field, typename GImpl>
class TCovariantLaplacian: public Module<CovariantLaplacianPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TCovariantLaplacian(const std::string name);
    // destructor
    virtual ~TCovariantLaplacian(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionCovariantLaplacian, 
                    ARG(TCovariantLaplacian<FIMPL::FermionField, GIMPL>), MAction);

/******************************************************************************
 *                      TCovariantLaplacian implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
TCovariantLaplacian<Field, GImpl>::TCovariantLaplacian(const std::string name)
: Module<CovariantLaplacianPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field, typename GImpl>
std::vector<std::string> TCovariantLaplacian<Field, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename Field, typename GImpl>
std::vector<std::string> TCovariantLaplacian<Field, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TCovariantLaplacian<Field, GImpl>::setup(void)
{
    LOG(Message) << "Setting up covariant Laplacian operator with gauge field '"
                 << par().gauge << "' and m^2= " << par().m2 << std::endl;
    auto &U = envGet(GaugeField, par().gauge);
    envCreateDerived(LinearOperatorBase<Field>, ARG(CovariantLaplacian<Field, GImpl>), 
                     getName(), 1, U, par().m2);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename GImpl>
void TCovariantLaplacian<Field, GImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_CovariantLaplacian_hpp_
