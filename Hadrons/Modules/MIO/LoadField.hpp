/*
 * LoadField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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
#ifndef Hadrons_MIO_LoadField_hpp_
#define Hadrons_MIO_LoadField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/FieldIo.hpp>
#include <Hadrons/Modules/MIO/SaveField.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Module to load a single field from disk                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadFieldPar,
                                    std::string, name,
                                    unsigned int, Ls,
                                    std::string, fileStem);
};

template <typename Field, typename FieldIo = Field>
class TLoadField: public Module<LoadFieldPar>
{
public:
    typedef FieldReader<Field, FieldIo> Reader;
public:
    // constructor
    TLoadField(const std::string name);
    // destructor
    virtual ~TLoadField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadPropagator, TLoadField<FIMPL::PropagatorField>, MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadPropagatorIo32, ARG(TLoadField<FIMPL::PropagatorField, FIMPLF::PropagatorField>), MIO);
#endif
MODULE_REGISTER_TMP(LoadColourMatrixField, TLoadField<GIMPL::GaugeLinkField>, MIO);

/******************************************************************************
 *                 TLoadField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
TLoadField<Field, FieldIo>::TLoadField(const std::string name)
: Module<LoadFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
std::vector<std::string> TLoadField<Field, FieldIo>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Field, typename FieldIo>
std::vector<std::string> TLoadField<Field, FieldIo>::getOutput(void)
{
    std::vector<std::string> out = {par().name};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
void TLoadField<Field, FieldIo>::setup(void)
{
    GridBase *grid = nullptr, *gridIo = nullptr;

    if (par().Ls > 1)
    {
        grid = envGetGrid(Field, par().Ls);
        if (!sameType<Field, FieldIo>())
        {
            gridIo = envGetGrid(FieldIo, par().Ls);
        }
        envCreateLat(Field, par().name, par().Ls);
    }
    else
    {
        grid = envGetGrid(Field);
        if (!sameType<Field, FieldIo>())
        {
            gridIo = envGetGrid(FieldIo);
        }
        envCreateLat(Field, par().name);
    }
    envTmp(Reader, "reader", 1, grid, gridIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
void TLoadField<Field, FieldIo>::execute(void)
{
    std::string          xmlString;
    auto                 &field = envGet(Field, par().name);
    GenericFieldMetadata md;
    envGetTmp(Reader, reader);

    LOG(Message) << "Loading field '" << par().name << "' from file '"
                 <<  resultFilename(par().fileStem, "bin") << "'" << std::endl;
    LOG(Message) << "Field type: " << typeName<Field>() << std::endl;
    if (!sameType<Field, FieldIo>())
    {
        LOG(Message) << "I/O type  : " << typeName<FieldIo>() << std::endl;
    }
    reader.open(resultFilename(par().fileStem, "bin"));
    reader.readHeader(xmlString);
    reader.readField(field, md);
    reader.close();
    LOG(Message) << "Field metadata:" << std::endl;
    LOG(Message) << "* name      : " << md.name << std::endl;
    LOG(Message) << "* moduleName: " << md.moduleName << std::endl;
    LOG(Message) << "* moduleType: " << md.moduleType << std::endl;
    LOG(Message) << "* moduleXml : " << md.moduleXml << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadField_hpp_
