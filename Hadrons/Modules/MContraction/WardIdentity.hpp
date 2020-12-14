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
                                    bool,        test_axial,
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
                                        std::vector<Scalar>,  DmuJmu,
                                        std::vector<Scalar>,  PDmuAmu,
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
    envTmpLat(PropagatorField, "vector_WI");
    envTmpLat(ComplexField, "vector_current");
    if (par().test_axial)
    {
        envTmpLat(PropagatorField, "psi");
        envTmpLat(ComplexField, "axial_defect");
        envTmpLat(ComplexField, "axial_defect_fin");
        envTmpLat(ComplexField, "PJ5q");
        envTmpLat(ComplexField, "PP");
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

    // Create results = zero
    Result result;
    result.mass = par().mass;
    const int nt { env().getDim(Tp) };
    result.DmuJmu.resize(nt, 0.);
    result.PDmuAmu.resize(nt, 0.);
    result.PP.resize(nt, 0.);
    result.PJ5q.resize(nt, 0.);
    result.mres.resize(nt, 0.);

    // Compute D_mu V_mu (D here is backward derivative)
    LOG(Message) << "Getting vector conserved current" << std::endl;
    envGetTmp(PropagatorField, tmp);
    envGetTmp(PropagatorField, vector_WI);
    vector_WI = Zero();
    for (unsigned int mu = 0; mu < Nd; ++mu)
    {
        act.ContractConservedCurrent(q, q, tmp, phys_source, Current::Vector, mu);
        tmp -= Cshift(tmp, mu, -1);
        vector_WI += tmp;
    }
    // Log ward identity D_mu V_mu = 0;
    LOG(Message) << "Vector Ward Identity check Delta_mu V_mu = " << norm2(vector_WI) << std::endl;

    // Save the spatial sum for each time-plane
    std::vector<TensorScalar> vector_buf;
    envGetTmp(ComplexField, vector_current);
    vector_current = trace(vector_WI);
    sliceSum(vector_current, vector_buf, Tp);
    for (int t = 0; t < nt; ++t)
    {
        result.DmuJmu[t] = TensorRemove(vector_buf[t]);
    }

    // Not sure why axial tests should be optional
    if (par().test_axial)
    {
        // Compute <P|D_mu A_mu>, D is backwards derivative.
        envGetTmp(PropagatorField, psi);
        envGetTmp(ComplexField, axial_defect);
        axial_defect = Zero();
        for (unsigned int mu = 0; mu < Nd; ++mu)
        {
            act.ContractConservedCurrent(q, q, tmp, phys_source, Current::Axial, mu);
            tmp -= Cshift(tmp, mu, -1);
            axial_defect += trace(g5 * tmp);
        }

        // Get <P|J5q> for 5D (zero for 4D) and <P|P>.
        envGetTmp(ComplexField, PP);
        envGetTmp(ComplexField, PJ5q);
        if (Ls_ > 1)
        {
            // <P|P>
            ExtractSlice(tmp, q, 0, 0);
            psi = 0.5 * (g5 * tmp - tmp);
            ExtractSlice(tmp, q, Ls_ - 1, 0);
            psi += 0.5 * (tmp + g5 * tmp);
            PP = trace(adj(psi) * psi);
            // <P|5Jq>
            act.ContractJ5q(q, PJ5q);
        }
        else
        {
            // 4d action
            PP = trace(adj(q) * q);
            PJ5q = Zero();
        }
        envGetTmp(ComplexField, axial_defect_fin);
        axial_defect_fin = axial_defect - 2. * (par().mass * PP + PJ5q);

        // Test ward identity <P|D_mu A_mu> = 2m<P|P> + 2<P|J5q>
        LOG(Message) << "|D_mu A_mu|^2 = " << norm2(axial_defect) << std::endl;
        LOG(Message) << "|PP|^2        = " << norm2(PP) << std::endl;
        LOG(Message) << "|PJ5q|^2      = " << norm2(PJ5q) << std::endl;
        LOG(Message) << "Axial Ward Identity defect Delta_mu A_mu = "
                     << norm2(axial_defect_fin) << std::endl;

        // Axial defect by timeslice.
        LOG(Message) << "Check Axial defect by timeslice" << std::endl;
        std::vector<TensorScalar> fin_buf;
        std::vector<TensorScalar> axial_buf;
        std::vector<TensorScalar> PP_buf;
        std::vector<TensorScalar> PJ5q_buf;
        sliceSum(axial_defect_fin, fin_buf, Tp);
        sliceSum(axial_defect, axial_buf, Tp);
        sliceSum(PP, PP_buf, Tp);
        sliceSum(PJ5q, PJ5q_buf, Tp);
        for (int t = 0; t < nt; ++t)
        {
            // This output can be compared with Grid Test_cayley_mres
            LOG(Message) << "t=" << t << ", Axial defect PAc=" << TensorRemove(axial_buf[t])
                         << ", PJ5q[" << t << "]=" << TensorRemove(fin_buf[t])
                         << ", PCAC_relation[" << t << "]=" << TensorRemove(fin_buf[t]) << std::endl;
            result.PDmuAmu[t] = TensorRemove(axial_buf[t]);
            result.PP[t]     = TensorRemove(PP_buf[t]);
            result.PJ5q[t]   = TensorRemove(PJ5q_buf[t]);
            result.mres[t]   = TensorRemove(PJ5q_buf[t]) / TensorRemove(PP_buf[t]);
        }
    }

    LOG(Message) << "Writing results to " << par().output << "." << std::endl;
    saveResult(par().output, "wardIdentity", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_WardIdentity_hpp_
