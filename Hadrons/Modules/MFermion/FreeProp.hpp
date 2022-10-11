/*
 * FreeProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
 * Author: guelpers <Vera.Guelpers@ed.ac.uk>
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


#ifndef Hadrons_MFermion_FreeProp_hpp_
#define Hadrons_MFermion_FreeProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Free fermion propagator                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class FreePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FreePropPar,
                                    std::string, source,
				                    std::string,  action,
				                    double, mass,
                                    std::string , boundary,
				                    std::string,  twist,
                                    std::string, outputTrace);
};

template <typename FImpl>
class TFreeProp: public Module<FreePropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename SpinMatrixField::scalar_object SpinMatrix;
    typedef Correlator<std::string, Complex>        ScalarResult;
    typedef Correlator<std::string, SpinMatrix>     SpinResult;
    class Result: Serializable
    {
        public:
            GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                            double,                    mass,
                                            std::vector<ComplexD>,     boundary,
                                            std::vector<double>,       twist,
                                            SpinResult,                full,
                                            std::vector<ScalarResult>, tr);
    };
public:
    // constructor
    TFreeProp(const std::string name);
    // destructor
    virtual ~TFreeProp(void) {};
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
};

MODULE_REGISTER_TMP(FreeProp, TFreeProp<FIMPL>, MFermion);

/******************************************************************************
 *                         TFreeProp implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TFreeProp<FImpl>::TFreeProp(const std::string name)
: Module<FreePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TFreeProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().action};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TFreeProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFreeProp<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().action);
    envCreateLat(PropagatorField, getName());
    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "source", Ls_);
        envTmpLat(FermionField, "sol", Ls_);
    }
    else
    {
       envTmpLat(FermionField, "source");
       envTmpLat(FermionField, "sol");
    }
    envTmpLat(FermionField, "tmp");
    if (Ls_ > 1)
    {
        envCreateLat(PropagatorField, getName() + "_5d", Ls_);
    } 
    if (!par().outputTrace.empty())
    {
        envTmpLat(ComplexField, "c");
        envTmpLat(SpinMatrixField, "cSpin");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFreeProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing free fermion propagator '" << getName() << "'"
                 << std::endl;
    
    std::string propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    auto        &prop    = envGet(PropagatorField, propName);
    auto        &fullSrc = envGet(PropagatorField, par().source);
    auto        &mat = envGet(FMat, par().action);
    RealD mass = par().mass;
    
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    LOG(Message) << "Calculating a free Propagator with mass " << mass 
		         << " using the action '" << par().action
                 << "' on source '" << par().source << "'" << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        LOG(Message) << "Calculation for spin= " << s << ", color= " << c
                     << std::endl;
        // source conversion for 4D sources
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
               PropToFerm<FImpl>(source, fullSrc, s, c);
            }
            else
            {
                PropToFerm<FImpl>(tmp, fullSrc, s, c);
                mat.ImportPhysicalFermionSource(tmp, source);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm<FImpl>(source, fullSrc, s, c);
            }
        }
        sol = Zero();
        std::vector<double> twist = strToVec<double>(par().twist);
        if(twist.size() != Nd)
        {
            HADRONS_ERROR(Size, "number of twist angles does not match number of dimensions");
        }
        std::vector<Complex> boundary = strToVec<Complex>(par().boundary);
        if(boundary.size() != Nd)
        {
            HADRONS_ERROR(Size, "number of boundary conditions does not match number of dimensions");
        }
	    mat.FreePropagator(source,sol,mass,boundary,twist);
        FermToProp<FImpl>(prop, sol, s, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            PropagatorField &p4d = envGet(PropagatorField, getName());
            mat.ExportPhysicalFermionSolution(sol, tmp);
            FermToProp<FImpl>(p4d, tmp, s, c);
        }
    }
    if (!par().outputTrace.empty())
    {
        PropagatorField       &p4d = envGet(PropagatorField, getName());
        Gamma                 gt(Gamma::Algebra::GammaT);
        Result                result;
        std::vector<TComplex> buf;

        envGetTmp(ComplexField, c);
        envGetTmp(SpinMatrixField, cSpin);
        cSpin = peekColour(p4d, 0, 0);
        result.mass     = par().mass;
        result.boundary = strToVec<ComplexD>(par().boundary);
        result.twist    = strToVec<double>(par().twist);
        sliceSum(cSpin, result.full.corr, Tp);
        result.full.info = "full spin";
        c = trace(p4d);
        sliceSum(c, buf, Tp);
        result.tr.resize(2);
        result.tr[0].info = "tr(S)";
        result.tr[0].corr.resize(buf.size());
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result.tr[0].corr[t] = TensorRemove(buf[t]);
        }
        c = trace(gt*p4d);
        sliceSum(c, buf, Tp);
        result.tr[1].info = "tr(GammaT*S)";
        result.tr[1].corr.resize(buf.size());
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result.tr[1].corr[t] = TensorRemove(buf[t]);
        }
        saveResult(par().outputTrace, "freeProp", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_FreeProp_hpp_
