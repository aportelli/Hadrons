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
                                        std::vector<Complex>, DmuPAmu, // D_mu trace(pseudoscalar*conserved_axial_mu)
                                        std::vector<Complex>, PP,      // Pseudoscalar density
                                        std::vector<Complex>, PJ5q,    // Midpoint axial current density
                                        std::vector<Complex>, mres,    // residual mass = PJ5q / PP
                                        std::vector<Complex>, VDmuJmu, // D_mu trace(local_vector*conserved_vector_mu)
                                        std::vector<Complex>, DefectPA); // DmuPAmu[t] - 2.*(result.mass*result.PP[t] + result.PJ5q[t])
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
    void SliceOut(std::vector<Complex> &Out, SlicedComplex &Sum, const ComplexField &f, bool bDiff=true) const
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
    // For 5d actions, I'll also need the 4d propagator, so we can compute pseudoscalar density
    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "ferm5d", Ls_); // One spin and colour of the 5d propagator
        envTmpLat(FermionField, "ferm4d"); // One spin and colour of the 4d propagator
        envTmpLat(PropagatorField, "psi"); // This will hold the 4d version of the 5d propagator
    }
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
    result.DmuPAmu.resize(nt, 0.);
    result.PP.resize(nt, 0.);
    result.PJ5q.resize(nt, 0.);
    result.mres.resize(nt, 0.);
    result.VDmuJmu.resize(nt, 0.);
    result.DefectPA.resize(nt, 0.);

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
    SliceOut(result.DmuJmu, sumSV, tmp_current);
    // Vector-vector current density
    tmp_current = trace(gT*tmp);
    SliceOut(result.VDmuJmu, sumVV, tmp_current);
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
        SliceOut(result.DmuPAmu, sumPA, tmp_current);
        // <P|J5q>
        act.ContractJ5q(prop, tmp_current);
        SlicedComplex sumPJ5q(nt);
        SliceOut(result.PJ5q, sumPJ5q, tmp_current, false);
        // <P|P>
        LOG(Message) << "Getting 4d propagator for " << par().prop << std::endl;
        envGetTmp(FermionField, ferm5d);
        envGetTmp(FermionField, ferm4d);
        envGetTmp(PropagatorField, psi);
        for (int s = 0; s < Ns; ++s)
        {
            for (int c = 0; c < Nc; ++c)
            {
                PropToFerm<FImpl>(ferm5d,prop,s,c);
                act.ExportPhysicalFermionSolution(ferm5d,ferm4d);
                FermToProp<FImpl>(psi,ferm4d,s,c);
            }
        }
        LOG(Message) << "Getting pseudoscalar density" << std::endl;
        tmp_current = trace(adj(psi) * psi);
        SlicedComplex sumPP(nt);
        SliceOut(result.PP, sumPP, tmp_current, false);
#ifdef  COMPARE_Test_Cayley_mres
        LOG(Message) << "Axial Ward Identity by timeslice" << std::endl;
        LOG(Message) << "Mass=" << result.mass << std::endl;
#endif
        for (int t = 0; t < nt; ++t)
        {
            result.DefectPA[t] = result.DmuPAmu[t] - 2.*(result.mass*result.PP[t] + result.PJ5q[t]);
            result.mres[t]     = result.PJ5q[t] / result.PP[t];
#ifdef  COMPARE_Test_Cayley_mres
            // This output can be compared with Grid Test_cayley_mres
            LOG(Message) << " t=" << t << ", DmuPAmu=" << real(result.DmuPAmu[t]) << ", PP=" << real(result.PP[t]) << ", PJ5q=" << real(result.PJ5q[t]) << ", PCAC/AWI defect=" << real(result.DefectPA[t]) << std::endl;
#endif
        }
    }

    LOG(Message) << "Writing results to " << par().output << "." << std::endl;
    saveResult(par().output, "wardIdentity", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WardIdentity_hpp_
