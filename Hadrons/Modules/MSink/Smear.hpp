/*
 * Smear.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
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

#ifndef Hadrons_MSink_Smear_hpp_
#define Hadrons_MSink_Smear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 Smear                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class SmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SmearPar,
                                    std::string, q,
                                    std::string, sink);
};

template <typename Field>
class TSmear: public Module<SmearPar>
{
public:
    typedef std::vector<typename Field::scalar_object> SlicedField;
    typedef std::function<SlicedField (const Field &)> SinkFn;
public:
    // constructor
    TSmear(const std::string name);
    // destructor
    virtual ~TSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Smear, TSmear<FIMPL::PropagatorField> , MSink);

/******************************************************************************
 *                          TSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TSmear<Field>::TSmear(const std::string name)
: Module<SmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TSmear<Field>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().sink};
    
    return in;
}

template <typename Field>
std::vector<std::string> TSmear<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TSmear<Field>::setup(void)
{
    envCreate(SlicedField, getName(), 1, env().getDim(Tp));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TSmear<Field>::execute(void)
{
    LOG(Message) << "Sink smearing field '" << par().q
                 << "' using sink function '" << par().sink << "'."
                 << std::endl;

    auto &sink = envGet(SinkFn, par().sink);
    auto &q    = envGet(Field, par().q);
    auto &out  = envGet(SlicedField, getName());
    
    out = sink(q);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_Smear_hpp_
