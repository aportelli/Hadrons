/*
 * SourceSink.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
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

#ifndef Hadrons_MSink_SourceSink_hpp_
#define Hadrons_MSink_SourceSink_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 SourceSink                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class SourceSinkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SourceSinkPar,
                                    std::string, source);
};

template <typename FImpl>
class TSourceSink: public Module<SourceSinkPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    // constructor
    TSourceSink(const std::string name);
    // destructor
    virtual ~TSourceSink(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SourceSink,       TSourceSink<FIMPL>,        MSink);

/******************************************************************************
 *                          TSourceSink implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSourceSink<FImpl>::TSourceSink(const std::string name)
: Module<SourceSinkPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSourceSink<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSourceSink<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSourceSink<FImpl>::setup(void)
{
    envCreate(SinkFn, getName(), 1, nullptr);
}


// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSourceSink<FImpl>::execute(void)
{
    LOG(Message) << "Setting up sink function with source '" << par().source << "' as the sink" << std::endl;

    PropagatorField &source  = envGet(PropagatorField, par().source);
    
    auto sink = [this, source](const PropagatorField &field)
    {
        SlicedPropagator res;
        PropagatorField tmp = source*field;
        
        sliceSum(tmp, res, Tp);
        
        return res;
    };
    envGet(SinkFn, getName()) = sink;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_SourceSink_hpp_
