/*
 * Utils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Matthew Black    <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MGradientFlow_Utils_hpp_
#define Hadrons_MGradientFlow_Utils_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MGradientFlow)

class GaugeResult : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeResult,
                                    std::vector<double>,    flowtime,
                                    std::vector<double>,    plaquette,
                                    std::vector<double>,    rectangle,
                                    std::vector<double>,    clover,
                                    std::vector<double>,    topocharge,
                                    std::vector<double>,    action,
                                    std::vector<ComplexD>,  polyakovX,
                                    std::vector<ComplexD>,  polyakovY,
                                    std::vector<ComplexD>,  polyakovZ,
                                    std::vector<ComplexD>,  polyakovT);
};

// additional action(s) /////////////////////////////////////////////////////
template <class GImpl>
class ZeuthenGaugeAction {
private:
    typename WilsonLoops<GImpl>::StapleAndRectStapleAllWorkspace workspace;
public:
    INHERIT_GIMPL_TYPES(GImpl);
    
    double beta;
    SymanzikGaugeAction<GImpl> SG;

    ZeuthenGaugeAction(double b): beta(b),SG(SymanzikGaugeAction<GImpl>(b)) {};

    virtual std::string action_name(){return "ZeuthenGaugeAction";}

    virtual double S(const GaugeField &U) {
        return SG.S(U);
    };

    virtual std::string LogParameters(){
        std::stringstream sstream;
        sstream << GridLogMessage << "["<<action_name() <<"] beta: " << beta << std::endl;
        return sstream.str();
    }


    virtual void deriv(const GaugeField &Umu, GaugeField &dSdU) {
                                                //  beta = 3.0, cl = -1.0/12.0 -> Symanzik
        double factor_p = 5.0/double(Nc)*0.5;   //   5.0 = beta*(1.0-8.0*cl)
        double factor_r = -0.25/double(Nc)*0.5; // -0.25 = beta*cl

        GridBase *grid = Umu.Grid();

        std::vector<GaugeLinkField> U (Nd,grid);

        for(int mu=0;mu<Nd;mu++){
            U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
        }

        std::vector<GaugeLinkField> RectStaple(Nd,grid), Staple(Nd,grid);
        WilsonLoops<GImpl>::StapleAndRectStapleAll(Staple, RectStaple, U, workspace);


        GaugeLinkField dSdU_mu(grid);
        GaugeLinkField tmq(grid),tmr(grid);

        for (int mu=0; mu < Nd; mu++){
            dSdU_mu = Ta(U[mu]*Staple[mu])*factor_p;
            dSdU_mu = dSdU_mu + Ta(U[mu]*RectStaple[mu])*factor_r;

            tmq = (adj(Cshift(U[mu],mu,-1)) * Cshift(dSdU_mu,mu,-1) * Cshift(U[mu],mu,-1));
            tmr = (U[mu] * Cshift(dSdU_mu,mu,1) * adj(U[mu]));

            dSdU_mu = 5.0/6.0*dSdU_mu + 1.0/12.0*tmq + 1.0/12.0*tmr;
            PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
        }
    };
};



// field evolution /////////////////////////////////////////////////////////////
template <typename FlowAction, typename GImpl, typename FImpl>
class Evolution {
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename GImpl::GaugeLinkField GaugeLinkField;

    private:
        GridBase *grid_;
        GaugeLinkField linkBuf_;
        GaugeField zBuf1_, zBuf2_, uBuf1_, uBuf2_;
        std::vector<GaugeField> Wi_;
        TimerArray &timer_;

    public:
        double epsilon, maxTau, taus;
        FlowAction SG;

        // constructor with beta
        Evolution(GridBase *grid, double beta, double step, double mTau, double ts, TimerArray *timer = nullptr) 
        : SG(FlowAction(beta)), epsilon(step), maxTau(mTau), taus(ts)
        , grid_(grid), linkBuf_(grid), zBuf1_(grid), zBuf2_(grid)
        , uBuf1_(grid), uBuf2_(grid), Wi_(5, grid), timer_(*timer)
        {};

        // constructor with c_plaq, c_rect
        Evolution(GridBase *grid, double c_plaq, double c_rect, double step, double mTau, double ts, TimerArray *timer = nullptr) 
        : SG(FlowAction(c_plaq, c_rect)), epsilon(step), maxTau(mTau), taus(ts)
        , grid_(grid), linkBuf_(grid), zBuf1_(grid), zBuf2_(grid)
        , uBuf1_(grid), uBuf2_(grid), Wi_(5, grid), timer_(*timer)
        {};

        // clover //////////////////////////////////////////////////////////////////////
        void siteClover(ComplexField &Clov, const GaugeField &U)
        {
            GaugeLinkField scaledUnit(U.Grid());
            Clov = Zero();
            scaledUnit = 1.0/Nc;
            for (int mu = 1; mu < Nd; mu++) {
                for (int nu = 0; nu < mu; nu++) {
                    linkBuf_ = PeekIndex<LorentzIndex>(U, mu);
                    WilsonLoops<GImpl>::FieldStrength(linkBuf_, U, mu, nu);
                    linkBuf_ -= trace(linkBuf_) * scaledUnit;
                    Clov -= trace(linkBuf_ * linkBuf_);
                }
            }
        }

        double avgClover(const GaugeField &Umu) 
        {
            ComplexField Clov(Umu.Grid());

            siteClover(Clov, Umu);
            auto Tc = sum(Clov);
            auto c = TensorRemove(Tc);

            double vol = Umu.Grid()->gSites();

            return c.real() / vol;
        }

        // polyakov loop in mu direction  //////////////////////////////////////////////
        ComplexD avgPolyakovLoopMu(const GaugeField &Umu, int mu) { // assuming Nd=4
            GaugeLinkField P(Umu.Grid());
            ComplexD out;

            double vol = Umu.Grid()->gSites();

            linkBuf_ = peekLorentz(Umu,mu);
            P = linkBuf_;
            for (int t=1; t < Umu.Grid()->GlobalDimensions()[mu]; t++) {
                P = GImpl::CovShiftForward(linkBuf_,mu,P);
            }
            RealD norm = 1.0/(Nc*vol);
            out = sum(trace(P))*norm;
            return out;
        }


        void gauge_RK(GaugeField &U) {
            Wi_[0] = U;                                     // W0
            timer_.startTimer("evolution_deriv");
            SG.deriv(U, zBuf1_);       
            timer_.stopTimer("evolution_deriv");                         
            zBuf1_ *= 0.25;                                 // Z0 = 1/4 * F(U)
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);   // U = W1 = exp(ep*Z0)*W0
            Wi_[1] = U;                                     // W1

            zBuf1_ *= -17.0/8.0;
            timer_.startTimer("evolution_deriv");
            SG.deriv(U,uBuf1_); 
            timer_.stopTimer("evolution_deriv");
            zBuf1_ += uBuf1_;                                // -17/32*Z0 + Z1
            zBuf1_ *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);    // U = W2 = exp(ep*Z)*W1
            Wi_[2] = U;                                      // W2

            zBuf1_ *= -4.0/3.0;
            timer_.startTimer("evolution_deriv");
            SG.deriv(U,uBuf1_); 
            timer_.stopTimer("evolution_deriv");
            zBuf1_ += uBuf1_;                                // 4/3*(17/36*Z0 -8/9*Z1) + Z2
            zBuf1_ *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
            Wi_[3] = U;                                      // W3
        };

        void gauge_RK_adaptive(GaugeField &U) {


            if (maxTau - taus < epsilon){
                epsilon = maxTau-taus;
            }
            uBuf2_ = U;
            Wi_[0] = U;                                        // W0
            SG.deriv(U, zBuf1_);
            zBuf2_ = -zBuf1_;
            zBuf1_ *= 0.25;                                    // Z0 = 1/4 * F(U)
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);      // U = W1 = exp(ep*Z0)*W0
            Wi_[1] = U;                                        // W1

            zBuf1_ *= -17.0/8.0;
            SG.deriv(U, uBuf1_); zBuf1_ += uBuf1_;             // -17/32*Z0 +Z1
            zBuf2_ += 2.0*uBuf1_;
            zBuf1_ *= 8.0/9.0;                                 // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);      // U = W2 = exp(ep*Z)*W1
            Wi_[2] = U;                                        // W2

            zBuf1_ *= -4.0/3.0;
            SG.deriv(U, uBuf1_); zBuf1_ += uBuf1_;             // 4/3*(17/36*Z0 -8/9*Z1) +Z2
            zBuf1_ *= 3.0/4.0;                                 // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(zBuf1_, U, -2.0*epsilon);      // V(t+e) = exp(ep*Z)*W2      
            Wi_[3] = U;                                        // W3
            
            GImpl::update_field(zBuf2_, uBuf2_, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
            Wi_[4] = uBuf2_;                                   // Uprime
        };

        void adaptive_eps(const GaugeField& U, const GaugeField& Uprime) {
            // Compute distance as norm^2 of the difference
            uBuf1_ = U - Uprime;
            double diff = norm2(uBuf1_);   
            // adjust integration step  

            taus += epsilon;
            epsilon = epsilon*0.95*std::pow(1e-4/diff,1./3.);
        };

        void evolve_gauge(GaugeField &U) {
            gauge_RK(U);
            U = Wi_[3];
        };
        
        void evolve_gauge_adaptive(GaugeField &U) {
            gauge_RK_adaptive(U);
            adaptive_eps(Wi_[3],Wi_[4]);
            U = Wi_[3];
        };

        
        void gauge_apply_boundary(GaugeField &Umu, std::vector<int> bc) {
            GaugeLinkField tmp1(Umu.Grid());
            GaugeLinkField tmp2(Umu.Grid());
            GaugeLinkField tmp3(Umu.Grid());
            Lattice<iScalar<vInteger>> coord(Umu.Grid());

            for (int mu = 0; mu < Nd; mu++) {
                LatticeCoordinate(coord,mu);
                
                tmp1 = PeekIndex<LorentzIndex>(Umu,mu);
                tmp2 = (double)bc[mu]*tmp1;
                int dimSize = Umu.Grid()->GlobalDimensions()[mu] - 1;
                tmp3 = where((coord == dimSize), tmp2, tmp1);
                PokeIndex<LorentzIndex>(Umu, tmp3, mu);
            }
        };

        PropagatorField generic_laplace(double a, double b, GaugeField &Umu, const PropagatorField& x_in, int skip_axis) {
            double Nx = Nd;
            if (skip_axis != -1) Nx--;

            PropagatorField x_out = (a + -2.0*Nx*b) * x_in;
            for (int mu = 0; mu < Nd; mu++) {
                if (mu != skip_axis) {
                    GaugeLinkField U = PeekIndex<LorentzIndex>(Umu, mu);
                    x_out += b*(GImpl::CovShiftForward(U,mu,x_in) + GImpl::CovShiftBackward(U,mu,x_in));
                }
            }
            return x_out;
        };

        void laplace_flow(GaugeField &W0, GaugeField &W1, GaugeField &W2, PropagatorField &prop) {
            PropagatorField psi1 = prop + (epsilon/4.0)*generic_laplace(0.0, 1.0, W0, prop, -1);
            PropagatorField psi2 = prop + (8.0*epsilon/9.0)*generic_laplace(0.0, 1.0, W1, psi1, -1) - (2.0*epsilon/9.0)*generic_laplace(0.0, 1.0, W0, prop, -1);
            PropagatorField psi3 = psi1 + (3.0*epsilon/4.0)*generic_laplace(0.0, 1.0, W2, psi2, -1);

            prop = psi3;
        };

        std::vector<GaugeField> & evolve_gaugeFF(GaugeField &U, std::vector<int> &bc) {
            gauge_RK(U);
            U = 1.0*Wi_[3];

            gauge_apply_boundary(Wi_[0],bc);
            gauge_apply_boundary(Wi_[1],bc);
            gauge_apply_boundary(Wi_[2],bc);

            return Wi_;
        };

        // gauge field status //////////////////////////////////////////////////////////
        void gauge_status(GaugeField &Umu, GaugeResult &result, double flowt)
        {
            double Q = WilsonLoops<GImpl>::TopologicalCharge(Umu);
            double plaq = WilsonLoops<GImpl>::avgPlaquette(Umu);
            double rect = WilsonLoops<GImpl>::avgRectangle(Umu);
            double clov = avgClover(Umu);
            double act = SG.S(Umu);
            ComplexD polyX = avgPolyakovLoopMu(Umu,0);
            ComplexD polyY = avgPolyakovLoopMu(Umu,1);
            ComplexD polyZ = avgPolyakovLoopMu(Umu,2);
            ComplexD polyT = avgPolyakovLoopMu(Umu,3);

            result.flowtime.push_back(flowt);
            result.plaquette.push_back(plaq);
            result.rectangle.push_back(rect);
            result.clover.push_back(clov);
            result.topocharge.push_back(Q);
            result.action.push_back(act);
            result.polyakovX.push_back(polyX);
            result.polyakovY.push_back(polyY);
            result.polyakovZ.push_back(polyZ);
            result.polyakovT.push_back(polyT);
        };
};

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_Utils_hpp_
