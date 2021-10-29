/*
 * SaveField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_SaveField_hpp_
#define Hadrons_MIO_SaveField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/FieldIo.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Module to save a single field to disk                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class GenericFieldMetadata: Serializable
{
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GenericFieldMetadata,
                                        std::string, name,
                                        std::string, moduleName,
                                        std::string, moduleType,
                                        std::string, moduleXml);
};

class SaveFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveFieldPar,
                                    std::string, name,
                                    std::string, fileStem);
};

template <typename Field, typename FieldIo = Field>
class TSaveField: public Module<SaveFieldPar>
{
public:
    typedef FieldWriter<Field, FieldIo> Writer;
public:
    // constructor
    TSaveField(const std::string name);
    // destructor
    virtual ~TSaveField(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SavePropagator, TSaveField<FIMPL::PropagatorField>, MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(SavePropagatorIo32, ARG(TSaveField<FIMPL::PropagatorField, FIMPLF::PropagatorField>), MIO);
#endif
MODULE_REGISTER_TMP(SaveColourMatrixField, TSaveField<GIMPL::GaugeLinkField>, MIO);

/******************************************************************************
 *                 TSaveField implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
TSaveField<Field, FieldIo>::TSaveField(const std::string name)
: Module<SaveFieldPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
std::vector<std::string> TSaveField<Field, FieldIo>::getInput(void)
{
    std::vector<std::string> in = {par().name};
    
    return in;
}

template <typename Field, typename FieldIo>
std::vector<std::string> TSaveField<Field, FieldIo>::getOutput(void)
{
    std::vector<std::string> out;
    
    return out;
}

template <typename Field, typename FieldIo>
std::vector<std::string> TSaveField<Field, FieldIo>::getOutputFiles(void)
{
    std::vector<std::string> out = {resultFilename(par().fileStem, "bin")};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
void TSaveField<Field, FieldIo>::setup(void)
{
    unsigned int Ls = env().getObjectLs(par().name);
    GridBase *grid = nullptr, *gridIo = nullptr;

    if (Ls > 1)
    {
        grid = envGetGrid(Field, Ls);
        if (!sameType<Field, FieldIo>())
        {
            gridIo = envGetGrid(FieldIo, Ls);
        }
    }
    else
    {
        grid = envGetGrid(Field);
        if (!sameType<Field, FieldIo>())
        {
            gridIo = envGetGrid(FieldIo);
        }
    }
    envTmp(Writer, "writer", 1, grid, gridIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field, typename FieldIo>
void TSaveField<Field, FieldIo>::execute(void)
{
    GenericFieldMetadata md;
    auto                 &field = envGet(Field, par().name);
    XmlWriter            xmlWriter("", "hadronsFieldIo");
    ModuleBase           *mod = vm().getModule(env().getObjectModule(par().name));
    envGetTmp(Writer, writer);

    LOG(Message) << "Saving field '" << par().name << "' to file '"
                 <<  resultFilename(par().fileStem, "bin") << "'" << std::endl;
    LOG(Message) << "Field type: " << typeName<Field>() << std::endl;
    if (!sameType<Field, FieldIo>())
    {
        LOG(Message) << "I/O type  : " << typeName<FieldIo>() << std::endl;
    }
    md.name       = par().name;
    md.moduleName = mod->getName();
    md.moduleType = vm().getModuleType(mod->getName());
    md.moduleXml  = mod->parString();
    write(xmlWriter, "originalFilename", resultFilename(par().fileStem, "bin"));
    write(xmlWriter, "ioModuleType", vm().getModuleType(getName()));
    writer.open(resultFilename(par().fileStem, "bin"));
    writer.writeHeader(xmlWriter, "hadronsFieldIo");
    writer.writeField(field, md);
    writer.close();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveField_hpp_
