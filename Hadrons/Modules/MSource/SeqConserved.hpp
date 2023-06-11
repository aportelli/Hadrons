/*
 * SeqConserved.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef Hadrons_MSource_SeqConserved_hpp_
#define Hadrons_MSource_SeqConserved_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source with insertion of conserved current. 
 Additionally optional insertion of a photon field A_\mu(x).
 -----------------------------
 * src_x = sum_{mu=mu_min}^{mu_max} 
     q_x * theta(x_3 - tA) * theta(tB - x_3) * J_mu * exp(i x.mom) (* A_\mu(x))
 
 * options:
 - q: input propagator (string)
 - action: fermion action used for propagator q (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - curr_type: type of conserved current to insert (Current)
 - mu_min: begin Lorentz Index (integer)
 - mu_max: end Lorentz Index (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 - photon: optional photon field (string)
 
 */

/******************************************************************************
 *                              SeqConserved                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SeqConservedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SeqConservedPar,
                                    std::string,  q,
                                    std::string,  action,
                                    std::string,  source,
                                    unsigned int, tA,
                                    unsigned int, tB,
                                    Current,      curr_type,
                                    unsigned int, mu_min,
                                    unsigned int, mu_max,
                                    std::string,  mom,
                                    std::string,  photon);
};

template <typename FImpl>
class TSeqConserved: public Module<SeqConservedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    typedef PhotonR::GaugeField     EmField;
public:
    // constructor
    TSeqConserved(const std::string name);
    // destructor
    virtual ~TSeqConserved(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void makeSource(PropagatorField &src, PropagatorField &q, PropagatorField &physSrc);
private:
    bool        SeqhasPhase_{false}; 
    std::string SeqmomphName_;
};

MODULE_REGISTER_TMP(SeqConserved, TSeqConserved<FIMPL>, MSource);
MODULE_REGISTER_TMP(ZSeqConserved, TSeqConserved<ZFIMPL>, MSource);


/******************************************************************************
 *                      TSeqConserved implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSeqConserved<FImpl>::TSeqConserved(const std::string name)
: Module<SeqConservedPar>(name)
, SeqmomphName_ (name + "_Seqmomph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSeqConserved<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().action, par().source};
    if (!par().photon.empty()) in.push_back(par().photon);
        
    return in;
}

template <typename FImpl>
std::vector<std::string> TSeqConserved<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
   return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqConserved<FImpl>::setup(void)
{
    auto Ls_ = env().getObjectLs(par().action);

    if (Ls_ > 1)
    {
        envTmpLat(PropagatorField, "src_tmp", Ls_);
    }
    else
    {
        envTmpLat(PropagatorField, "src_tmp");
    }
    if (envHasType(PropagatorField, par().q))
    {
        if (Ls_ > 1)
        {
            envCreateLat(PropagatorField, getName(), Ls_);
        }
        else
        {
            envCreateLat(PropagatorField, getName());
        }
    }
    else if (envHasType(std::vector<PropagatorField>, par().q))
    {
        auto &q = envGet(std::vector<PropagatorField>, par().q);

        if (Ls_ > 1)
        {
            envCreate(std::vector<PropagatorField>, getName(), Ls_, q.size(),
                      envGetGrid(PropagatorField, Ls_));
        }
        else
        {
            envCreate(std::vector<PropagatorField>, getName(), 1, q.size(),
                      envGetGrid(PropagatorField));
        }
    }
    else
    {
        HADRONS_ERROR_REF(ObjectType, "object '" + par().q 
                          + "' has an incompatible type ("
                          + env().getObjectType(par().q)
                          + ")", env().getObjectAddress(par().q))
    }
    envCacheLat(LatticeComplex, SeqmomphName_);
    envTmpLat(LatticeComplex, "coor");
    envTmpLat(LatticeComplex, "latt_compl");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSeqConserved<FImpl>::makeSource(PropagatorField &src, PropagatorField &q, PropagatorField &physSrc)
{
    auto &mat = envGet(FMat, par().action);

    envGetTmp(PropagatorField, src_tmp);
    envGetTmp(LatticeComplex, latt_compl);

    src     = Zero();
    src_tmp = src;
    //exp(ipx)
    auto &mom_phase = envGet(LatticeComplex, SeqmomphName_);
    if (!SeqhasPhase_)
    {    
        std::vector<Real> mom = strToVec<Real>(par().mom);
        mom_phase = Zero();
        Complex           i(0.0,1.0);
        envGetTmp(LatticeComplex, coor);
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            mom_phase = mom_phase + (mom[mu]/env().getDim(mu))*coor;
        }
        mom_phase = exp((Real)(2*M_PI)*i*mom_phase);
        SeqhasPhase_ = true;
    }
    LOG(Message) << "Inserting momentum " << strToVec<Real>(par().mom) << std::endl;
    if (!par().photon.empty())    	
    {
	 LOG(Message) << "Inserting the stochastic photon field " << par().photon << std::endl;
    }
    for(unsigned int mu=par().mu_min;mu<=par().mu_max;mu++)
    {
        if (!par().photon.empty())    	
        {
	    //Get the stochastic photon field, if required
            auto &stoch_photon = envGet(EmField,  par().photon);
    	    latt_compl =  PeekIndex<LorentzIndex>(stoch_photon, mu) * mom_phase;
        }
        else
        {
            latt_compl = mom_phase;
        } 

    	mat.SeqConservedCurrent(q, src_tmp, physSrc, par().curr_type, mu, 
                             par().tA, par().tB, latt_compl);
	    src += src_tmp;
    }	
}

template <typename FImpl>
void TSeqConserved<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating sequential source(s) with conserved "
                     << par().curr_type << " current for the action '"
                     << par().action << "'\nat " 
		             << "t = " << par().tA << " summed over the indices " 
		             << par().mu_min << " <= mu <= " << par().mu_max 
		             << std::endl;
    }
    else
    {
        LOG(Message) << "Generating sequential source(s) with conserved "
                     << par().curr_type << " current for the action '"
                     << par().action << "'\nfor " 
                     << par().tA << " <= t <= " 
                     << par().tB << " summed over the indices " 
		             << par().mu_min << " <= mu <= " << par().mu_max
	                 << std::endl;
    }

    if (envHasType(PropagatorField, par().q))
    {
        auto  &src     = envGet(PropagatorField, getName());
        auto  &physSrc = envGet(PropagatorField, par().source);
        auto  &q       = envGet(PropagatorField, par().q);

        LOG(Message) << "Using propagator '" << par().q << "'" << std::endl;
        makeSource(src, q, physSrc);
    }
    else
    {
        auto  &src     = envGet(std::vector<PropagatorField>, getName());
        auto  &physSrc = envGet(std::vector<PropagatorField>, par().source);
        auto  &q       = envGet(std::vector<PropagatorField>, par().q);

        for (unsigned int i = 0; i < q.size(); ++i)
        {
            LOG(Message) << "Using element " << i << " of propagator vector '" 
                         << par().q << "'" << std::endl;
            makeSource(src[i], q[i], physSrc[i]);
        }
    }
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SeqConserved_hpp_
