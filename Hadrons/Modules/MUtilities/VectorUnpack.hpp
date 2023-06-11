/*
 * VectorUnpack.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Michael Marshall <michael.marshall@ed.ac.uk>
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
#ifndef Hadrons_MUtilities_VectorUnpack_hpp_
#define Hadrons_MUtilities_VectorUnpack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                Utility module to unpack a vector of fields                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class VectorUnpackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VectorUnpackPar,
                                    std::string,  input,
                                    unsigned int, size,
                                    std::vector<std::string>, fields);
};

template <typename Field>
class TVectorUnpack: public Module<VectorUnpackPar>
{
public:
    // constructor
    TVectorUnpack(const std::string name);
    // destructor
    virtual ~TVectorUnpack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
protected:
    inline std::string getPackedName(unsigned int n) const
    {
        if (n < par().fields.size())
        {
            return par().fields[n];
        }
        return getName() + "_" + std::to_string(n);
    }
};

MODULE_REGISTER_TMP(ComplexVectorUnpack, TVectorUnpack<FIMPL::ComplexField>, MUtilities);
MODULE_REGISTER_TMP(FermionVectorUnpack, TVectorUnpack<FIMPL::FermionField>, MUtilities);
MODULE_REGISTER_TMP(PropagatorVectorUnpack, TVectorUnpack<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                       TVectorUnpack implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TVectorUnpack<Field>::TVectorUnpack(const std::string name)
: Module<VectorUnpackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TVectorUnpack<Field>::getInput(void)
{
    std::vector<std::string> in = {par().input};
    
    return in;
}

template <typename Field>
std::vector<std::string> TVectorUnpack<Field>::getOutput(void)
{
    std::vector<std::string> out;
    for (unsigned int i = 0; i < par().size; ++i)
    {
        out.push_back(getPackedName(i));
    }
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorUnpack<Field>::setup(void)
{
    auto         &vec  = envGet(std::vector<Field>, par().input);
    unsigned int Ls    = env().getObjectLs(par().input);
    auto         *grid = vec[0].Grid();

    if (vec.size() != par().size)
    {
        HADRONS_ERROR(Size,"Mismatch between vector size ("
                            + std::to_string(vec.size())
                            + ") and module parameter size ("
                            + std::to_string(par().size) + ").");
    }

    if (par().fields.size() && par().fields.size() != vec.size())
    {
        HADRONS_ERROR(Size,"Mismatch between vector size ("
                            + std::to_string(vec.size())
                            + ") and number of field names ("
                            + std::to_string(par().fields.size()) + ").");
    }

    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        envCreate(Field, getPackedName(i), Ls, grid);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorUnpack<Field>::execute(void)
{
    auto &vec = envGet(std::vector<Field>, par().input);

    LOG(Message) << "Unpacking vector '" << par().input << "'" << std::endl;
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
        auto &veci = envGet(Field, getPackedName(i));

        veci = vec[i];
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_VectorUnpack_hpp_
