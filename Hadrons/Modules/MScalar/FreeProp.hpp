/*
 * FreeProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MScalar_FreeProp_hpp_
#define Hadrons_MScalar_FreeProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MAction/Laplacian.hpp>
#include <Hadrons/Modules/MScalar/Scalar.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                               FreeProp                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalar)

class FreePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FreePropPar,
                                    std::string, source,
                                    double,      mass,
                                    std::string, output,
                                    bool, useFft,
                                    bool, cache);
};

template <typename SImpl>
class TFreeProp: public Module<FreePropPar>
{
public:
    typedef typename SImpl::Field                  Field;
    typedef MAction::Laplacian<Field>              LapMat;
    typedef HermitianLinearOperator<LapMat, Field> LapOp;
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
    std::string freeMomPropName_;
    bool        freePropDone_, freeMomPropDone_;
};

MODULE_REGISTER_TMP(FreeProp, TFreeProp<SIMPL>, MScalar);

/******************************************************************************
*                        TFreeProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TFreeProp<SImpl>::TFreeProp(const std::string name)
: Module<FreePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TFreeProp<SImpl>::getInput(void)
{
    return {par().source};
}

template <typename SImpl>
std::vector<std::string> TFreeProp<SImpl>::getOutput(void)
{
    return {getName(), getName()+"_sliceSum"};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TFreeProp<SImpl>::setup(void)
{
    if (par().useFft)
    {
        freeMomPropName_ = FREEMOMPROP(par().mass);
        freeMomPropDone_ = env().hasCreatedObject(freeMomPropName_);
        envCacheLat(Field, freeMomPropName_);
    }
    if (par().cache)
    {
        freePropDone_ = env().hasCreatedObject(getName());
        envCacheLat(Field, getName());
    }
    else
    {
        freePropDone_ = false;
        envCreateLat(Field, getName());
    }
    envCreate(HadronsSerializable, getName() + "_sliceSum", 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TFreeProp<SImpl>::execute(void)
{
    auto &prop   = envGet(Field, getName());
    auto &source = envGet(Field, par().source);

    if (!((par().cache) && freePropDone_))
    {
        LOG(Message) << "Computing free scalar propagator..." << std::endl;
        if (par().useFft)
        {
            auto &freeMomProp = envGet(Field, freeMomPropName_);
            if (!freeMomPropDone_)
            {
                LOG(Message) << "Caching momentum space free scalar propagator"
                            << " (mass= " << par().mass << ")..." << std::endl;
                SImpl::MomentumSpacePropagator(freeMomProp, par().mass);
            }
            SImpl::FreePropagator(source, prop, freeMomProp);
        }
        else
        {
            LapMat lap(par().mass*par().mass, getGrid<Field>());
            LapOp op(lap);
            ConjugateGradient<Field> cg(1.0e-8, 10000);

            cg(op, source, prop);
        }
    }
    else
    {
        LOG(Message) << "Using cached free scalar propagator" << std::endl;
    }
    
    std::vector<TComplex> buf;
    std::vector<Complex>  result;

    sliceSum(prop, buf, Tp);
    result.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result[t] = TensorRemove(buf[t]);
    }
    envGet(HadronsSerializable, getName()+"_sliceSum") = result;
    saveResult(par().output, "freeprop", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalar_FreeProp_hpp_
