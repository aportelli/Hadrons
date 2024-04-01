/*
 * Laplacian.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MAction_Laplacian_hpp_
#define Hadrons_MAction_Laplacian_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Laplacian operator                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class LaplacianPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LaplacianPar,
                                    double, m2);
};

template<class Field> 
class Laplacian : public SparseMatrixBase<Field>
{
public:
    typedef typename Field::vector_object vobj;
    typedef CartesianStencil<vobj, vobj, SimpleStencilParams> Stencil;
public:
    Laplacian(const double m2, GridBase *g)
    : m2_(m2), grid_(g)
    {
        auto nd = grid_->Nd();
        for (unsigned int mu = 0; mu < nd; ++mu)
        {
            dir_.push_back(mu);
            disp_.push_back(1);
        }
        for (unsigned int mu = 0; mu < nd; ++mu)
        {
            dir_.push_back(mu);
            disp_.push_back(-1);
        }
        st_.reset(new Stencil(grid_, 2*nd, Even, dir_, disp_, SimpleStencilParams()));
    };

    virtual GridBase *Grid(void) { return grid_; };
    virtual void M(const Field &_in, Field &_out)
    {
        auto nd = grid_->Nd();

        st_->HaloExchange(_in, comp_);

        auto st = st_->View(AcceleratorRead);
        auto buf = st_->CommBuf();
        autoView(in, _in, AcceleratorRead);
        autoView(out, _out, AcceleratorWrite);
        const int Nsimd = vobj::Nsimd();
        const uint64_t NN = grid_->oSites();

        typedef decltype(coalescedRead(in[0])) calcObj;

        accelerator_for(ss, NN, Nsimd,
        {
            StencilEntry *SE;
            const int lane = acceleratorSIMTlane(Nsimd);
            calcObj chi;
            calcObj res;
            int ptype;

            res = coalescedRead(in[ss]) * (m2_ + 2.0 * nd);
            for (unsigned int mu = 0; mu < 2 * nd; ++mu)
            {
                SE = st.GetEntry(ptype, mu, ss);
                if (SE->_is_local)
                {
                    int perm = SE->_permute;
                    chi = coalescedReadPermute(in[SE->_offset], ptype, perm, lane);
                }
                else
                {
                    chi = coalescedRead(buf[SE->_offset], lane);
                }
                acceleratorSynchronise();
                res -= chi;
            }
            coalescedWrite(out[ss], res, lane);
        });
    };
    virtual void Mdag(const Field &in, Field &out) { M(in, out); };
    virtual void Mdiag(const Field &in, Field &out) { assert(0); };
    virtual void Mdir(const Field &in, Field &out, int dir, int disp) { assert(0); }; 
    virtual void MdirAll(const Field &in, std::vector<Field> &out) { assert(0); };
private:
    double m2_;
    GridBase *grid_;
    std::unique_ptr<Stencil> st_;
    SimpleCompressor<vobj> comp_;
    std::vector<int> dir_, disp_;
};


template <typename Field>
class TLaplacian: public Module<LaplacianPar>
{
public:
    // constructor
    TLaplacian(const std::string name);
    // destructor
    virtual ~TLaplacian(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionLaplacian, TLaplacian<FIMPL::FermionField>, MAction);

/******************************************************************************
 *                      TLaplacian implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TLaplacian<Field>::TLaplacian(const std::string name)
: Module<LaplacianPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TLaplacian<Field>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field>
std::vector<std::string> TLaplacian<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TLaplacian<Field>::setup(void)
{
    LOG(Message) << "Setting up Laplacian operator with m^2= " << par().m2 << std::endl;
    envCreateDerived(SparseMatrixBase<Field>, ARG(Laplacian<Field>), 
                     getName(), 1, par().m2, getGrid<Field>());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TLaplacian<Field>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_Laplacian_hpp_
