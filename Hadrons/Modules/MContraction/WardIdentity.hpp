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
 - q:          propagator NB: "_5d" suffix will be appended (string)
 - source:     source module for the quark, used to remove contact terms (string)
 - action:     action module used for propagator solution (string)
 - mass:       mass of quark (double)
 - output:     filename for output (string)
 - test_axial: whether or not to test PCAC relation.
*/

/******************************************************************************
 *                              WardIdentity                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class WardIdentityPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WardIdentityPar,
                                    std::string, q,
                                    std::string, action,
                                    std::string, source,
                                    double,      mass,
                                    //bool,        test_axial,
                                    std::string, output);
};

template <typename FImpl>
class TWardIdentity: public Module<WardIdentityPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    using Scalar = typename FImpl::Scalar; // Effectively a Complex
    using TensorScalar = iScalar<iScalar<iScalar<Scalar>>>;
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        double,               mass,
                                        std::vector<Scalar>,  SV,
                                        std::vector<Scalar>,  VV,
                                        std::vector<Scalar>,  PA,
                                        std::vector<Scalar>,  DmuPA,
                                        std::vector<Scalar>,  DefectPA,
                                        std::vector<Scalar>,  PP,
                                        std::vector<Scalar>,  PJ5q,
                                        std::vector<Scalar>,  mres);
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
    void SliceOut(std::vector<Scalar> &Out, std::vector<TensorScalar> &Sum, ComplexField &f)
    {
        sliceSum(f, Sum, Tp);
        const auto nt = Sum.size();
        for (size_t t = 0; t < nt; ++t)
        {
            Out[t] = TensorRemove(Sum[t] - Sum[(t-1+nt)%nt]);
        }
    }
private:
    unsigned int Ls_;
    std::string qName;
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
    qName = par().q; // If you want the 5d version, you must say so
    return {qName, par().action, par().source};
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
    // The quark can be 4d or 5d, but must match the action
    const unsigned int ActionLs_{ env().getObjectLs(par().action) };
    Ls_ = env().getObjectLs( qName );
    if (Ls_ != ActionLs_)
    {
        std::string sError{ "Ls mismatch: quark Ls="};
        sError.append( std::to_string( Ls_ ) );
        sError.append( ", action Ls=" );
        sError.append( std::to_string( ActionLs_ ) );
        sError.append( ". Did you mean quark='" );
        sError.append( qName );
        sError.append( "_5d' ?" );
        HADRONS_ERROR(Size, sError);
    }
    // These temporaries are always 4d
    envTmpLat(PropagatorField, "tmp");
    envTmpLat(ComplexField, "tmp_current");
    //if (par().test_axial)
    {
        envTmpLat(PropagatorField, "psi");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWardIdentity<FImpl>::execute(void)
{
    LOG(Message) << "Performing Ward Identity checks for quark " << qName << std::endl;
    auto &q = envGet(PropagatorField, qName);
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
    result.SV.resize(nt, 0.);
    result.VV.resize(nt, 0.);
    result.PA.resize(nt, 0.);
    result.DmuPA.resize(nt, 0.);
    result.DefectPA.resize(nt, 0.);
    result.PP.resize(nt, 0.);
    result.PJ5q.resize(nt, 0.);
    result.mres.resize(nt, 0.);

    // Compute D_mu V_mu (D here is backward derivative)
    // There is no point performing Dmu on spatial directions, because after the spatial sum, these become zero
    envGetTmp(PropagatorField, tmp);
    envGetTmp(ComplexField, tmp_current);
    LOG(Message) << "Getting vector conserved current" << std::endl;
    act.ContractConservedCurrent(q, q, tmp, phys_source, Current::Vector, Tdir);
    // Scalar-vector current density
    tmp_current = trace(tmp);
    std::vector<TensorScalar> tmp_buf(nt);
    SliceOut(result.SV, tmp_buf, tmp_current);
    // Vector-vector current density
    tmp_current = trace(gT*tmp);
    SliceOut(result.VV, tmp_buf, tmp_current);
    // For comparison with Grid Test_Cayley_mres
    LOG(Message) << "Vector Ward Identity by timeslice" << std::endl;
    for (int t = 0; t < nt; ++t)
    {
        LOG(Message) << " t=" << t << ", SV=" << real(result.SV[t]) << ", VV=" << real(result.VV[t]) << std::endl;
    }

    // Save the spatial sum for each time-plane

    // Not sure why axial tests should be optional
    //if (par().test_axial)
    {
        LOG(Message) << "Getting axial conserved current" << std::endl;
        act.ContractConservedCurrent(q, q, tmp, phys_source, Current::Axial, Tdir);
        // Pseudoscalar-Axial current density
        tmp_current = trace(g5 * tmp);
        SliceOut(result.PA, tmp_buf, tmp_current);

        // Get <P|J5q> for 5D (zero for 4D) and <P|P>.
        if (Ls_ > 1)
        {
            // <P|5Jq>
            act.ContractJ5q(q, tmp_current);
            SliceOut(result.PJ5q, tmp_buf, tmp_current);
            // <P|P>
            envGetTmp(PropagatorField, psi);
            ExtractSlice(tmp, q, 0, 0);
            psi = 0.5 * (g5 * tmp - tmp);
            ExtractSlice(tmp, q, Ls_ - 1, 0);
            psi += 0.5 * (tmp + g5 * tmp);
            tmp_current = trace(adj(psi) * psi);
        }
        else
        {
            // 4d action
            tmp_current = trace(adj(q) * q);
        }
        SliceOut(result.PP, tmp_buf, tmp_current);
        LOG(Message) << "Axial Ward Identity by timeslice" << std::endl;
        for (int t = 0; t < nt; ++t)
        {
            result.DmuPA[t]    = result.PA[t] - result.PA[(t-1+nt)%nt];
            result.DefectPA[t] = result.DmuPA[t] - 2.*(result.mass*result.PP[t] + result.PJ5q[t]);
            result.mres[t]     = result.PJ5q[t] / result.PP[t];
            // This output can be compared with Grid Test_cayley_mres
            LOG(Message) << " t=" << t << ", PAc=" << real(result.PA[t]) << ", PJ5q=" << real(result.PJ5q[t])
                         << ", PCAC/AWI defect=" << real(result.DefectPA[t]) << std::endl;
        }
    }

    LOG(Message) << "Writing results to " << par().output << "." << std::endl;
    saveResult(par().output, "wardIdentity", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WardIdentity_hpp_
