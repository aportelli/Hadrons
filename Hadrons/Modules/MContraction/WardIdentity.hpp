/*
 * WardIdentity.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
 * Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
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
#ifndef Hadrons_MContraction_WardIdentity_hpp_
#define Hadrons_MContraction_WardIdentity_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
  Ward Identity contractions
 -----------------------------
 
 * options:
 - prop:       propagator. Must match the action, i.e. 5D action needs 5D propagator
 - action:     action module used for propagator solution (string)
 - source:     source module for the quark, used to remove contact terms (string)
 - mass:       mass of quark (double)
 - output:     filename for output (string)
*/

/******************************************************************************
 *                              WardIdentity                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class WardIdentityPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WardIdentityPar,
                                    std::string, prop,   // Name of the propagator we are checking Ward identity
                                    std::string, action,
                                    std::string, source,
                                    double,      mass,
                                    std::string, output);
};

template <typename FImpl>
class TWardIdentity: public Module<WardIdentityPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        double,               mass,
                                        std::vector<Complex>, DmuJmu,  // D_mu trace(Scalar*conserved_vector_mu)
                                        std::vector<Complex>, VDmuJmu, // D_mu trace(local_vector*conserved_vector_mu)
                                        std::vector<Complex>, PJ5q,    // Midpoint axial current density
                                        std::vector<Complex>, PA0);    //      trace(pseudoscalar*conserved_axial_0)
                                     // std::vector<Complex>, DmuPAmu, // D_mu trace(pseudoscalar*conserved_axial_mu)
                                     // std::vector<Complex>, PP,      // Pseudoscalar density
                                     // std::vector<Complex>, mres,    // residual mass = <PJ5q> / <PP>
                                     // std::vector<Complex>, DefectPA,// DmuPAmu[t] - 2.*(mass*PP[t] + PJ5q[t])
                                     // std::vector<Complex>, PALocal0,//trace( temporal local axial - pseudoscalar )
    };

public:
    // constructor
    TWardIdentity(const std::string name);
    // destructor
    virtual ~TWardIdentity(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // Perform Slice Sum and then save delta
    void SliceOut(std::vector<Complex> &Out, SlicedComplex &Sum, const ComplexField &f, bool bDiff) const
    {
        sliceSum(f, Sum, Tp);
        const auto nt = Sum.size();
        for (size_t t = 0; t < nt; ++t)
        {
            Out[t] = TensorRemove(bDiff ? Sum[t] - Sum[(t-1+nt)%nt] : Sum[t]);
        }
    }
private:
    unsigned int Ls_;
};

MODULE_REGISTER_TMP(WardIdentity, TWardIdentity<FIMPL>, MContraction);
MODULE_REGISTER_TMP(ZWardIdentity, TWardIdentity<ZFIMPL>, MContraction);

/******************************************************************************
 *                     TWardIdentity implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWardIdentity<FImpl>::TWardIdentity(const std::string name)
: Module<WardIdentityPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWardIdentity<FImpl>::getInput(void)
{
    return { par().prop, par().action, par().source };
}

template <typename FImpl>
std::vector<std::string> TWardIdentity<FImpl>::getOutput(void)
{
  return {};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWardIdentity<FImpl>::setup(void)
{
    // The propagator can be 4d or 5d, but must match the action
    const unsigned int ActionLs_{ env().getObjectLs(par().action) };
    Ls_ = env().getObjectLs( par().prop );
    if (Ls_ != ActionLs_)
    {
        std::string sError{ "Ls mismatch: propagator Ls="};
        sError.append( std::to_string( Ls_ ) );
        sError.append( ", action Ls=" );
        sError.append( std::to_string( ActionLs_ ) );
        HADRONS_ERROR(Size, sError);
    }
    // These temporaries are always 4d
    envTmpLat(PropagatorField, "tmp");
    envTmpLat(ComplexField, "tmp_current");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWardIdentity<FImpl>::execute(void)
{
    LOG(Message) << "Performing Ward Identity checks for propagator " << par().prop << std::endl;
    auto &prop = envGet(PropagatorField, par().prop);
    LOG(Message) << "Action " << par().action << std::endl;
    auto &act = envGet(FMat, par().action);
    LOG(Message) << "Physical source " << par().source << std::endl;
    auto &phys_source = envGet(PropagatorField, par().source);
    Gamma g5(Gamma::Algebra::Gamma5);
    Gamma gT(Gamma::Algebra::GammaT);

    // Create results = zero
    Result result;
    result.mass = par().mass;
    const int nt { env().getDim(Tp) };
    result.DmuJmu.resize(nt, 0.);
    result.VDmuJmu.resize(nt, 0.);
    result.PJ5q.resize(nt, 0.);
    result.PA0.resize(nt, 0.);

    // Compute D_mu V_mu (D here is backward derivative)
    // There is no point performing Dmu on spatial directions, because after the spatial sum, these become zero
    envGetTmp(PropagatorField, tmp);
    envGetTmp(ComplexField, tmp_current);
    SlicedComplex sumSV(nt);
    SlicedComplex sumVV(nt);
    LOG(Message) << "Getting vector conserved current" << std::endl;
    act.ContractConservedCurrent(prop, prop, tmp, phys_source, Current::Vector, Tdir);
    // Scalar-vector current density
    tmp_current = trace(tmp);
    SliceOut(result.DmuJmu, sumSV, tmp_current, true);
    // Vector-vector current density
    tmp_current = trace(gT*tmp);
    SliceOut(result.VDmuJmu, sumVV, tmp_current, true);
//#define COMPARE_Test_Cayley_mres
#ifdef  COMPARE_Test_Cayley_mres
    // For comparison with Grid Test_Cayley_mres
    LOG(Message) << "Vector Ward Identity by timeslice" << std::endl;
    for (int t = 0; t < nt; ++t)
    {
        LOG(Message) << " t=" << t << ", SV=" << real(TensorRemove(sumSV[t]))
                     << ", VV=" << real(TensorRemove(sumVV[t])) << std::endl;
    }
#endif

    // Test axial Ward identity for 5D actions
    if (Ls_ > 1)
    {
        LOG(Message) << "Getting axial conserved current" << std::endl;
        act.ContractConservedCurrent(prop, prop, tmp, phys_source, Current::Axial, Tdir);
        // Pseudoscalar-Axial current density
        tmp_current = trace(g5 * tmp);
        SlicedComplex sumPA(nt);
        // Save temporal component of pseudoscalar-(partially) conserved axial
        // \mathcal{A}_0 from eq (37) in https://arxiv.org/pdf/hep-lat/0612005.pdf
        SliceOut(result.PA0, sumPA, tmp_current, false);
        // <P|J5q>
        act.ContractJ5q(prop, tmp_current);
        SlicedComplex sumPJ5q(nt);
        SliceOut(result.PJ5q, sumPJ5q, tmp_current, false);
    }

    LOG(Message) << "Writing results to " << par().output << "." << std::endl;
    saveResult(par().output, "wardIdentity", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WardIdentity_hpp_
