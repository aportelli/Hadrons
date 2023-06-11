/*
 * StochasticQedL.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGauge_StochasticQedL_hpp_
#define Hadrons_MGauge_StochasticQedL_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EmField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Stochastic  QED_L                                        *
 *  Potentially with IR improvement coefficients                              *
 *  cf. https://journals.aps.org/prd/pdf/10.1103/PhysRevD.99.034510 (Eq. 85)  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGauge)

class StochasticQedLPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochasticQedLPar,
                                    QedGauge, gauge,
                                    std::string, improvement);
};

template <typename VType>
class TStochasticQedL: public Module<StochasticQedLPar>
{
public:
    typedef TEmFieldGenerator<VType>    EmGen;
    typedef typename EmGen::GaugeField  GaugeField;
    typedef typename EmGen::ScalarField ScalarField;
public:
    // constructor
    TStochasticQedL(const std::string name);
    // destructor
    virtual ~TStochasticQedL(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool weightDone_;
};

MODULE_REGISTER_TMP(StochasticQedL, ARG(TStochasticQedL<vComplex>), MGauge);

/******************************************************************************
 *                     TStochasticQedL implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename VType>
TStochasticQedL<VType>::TStochasticQedL(const std::string name)
: Module<StochasticQedLPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename VType>
std::vector<std::string> TStochasticQedL<VType>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename VType>
std::vector<std::string> TStochasticQedL<VType>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename VType>
void TStochasticQedL<VType>::setup(void)
{
    weightDone_ = env().hasCreatedObject("_" + getName() + "_weight");
    envCacheLat(ScalarField, "_" + getName() + "_weight");
    envCreateLat(GaugeField, getName());
    envTmp(EmGen, "gen", 1, envGetGrid(GaugeField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename VType>
void TStochasticQedL<VType>::execute(void)
{
    auto &a = envGet(GaugeField, getName());
    auto &w = envGet(ScalarField, "_" + getName() + "_weight");
    std::vector<double> improvement = strToVec<double>(par().improvement);
    envGetTmp(EmGen, gen);
    if (!weightDone_)
    {
        LOG(Message) << "Caching stochastic QED_L EM potential weights" << std::endl;
        if (improvement.size() > 0)
        {
            LOG(Message) << "Improvement coefficients " << improvement << std::endl;
        }
        gen.makeWeightsQedL(w, improvement);
    }
    LOG(Message) << "Generating stochastic EM potential (gauge: " << par().gauge << ")" << std::endl;
    auto tr = gen.getGaugeTranform(par().gauge);
    gen(a, rng4d(), w, tr);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGauge_StochasticQedL_hpp_
