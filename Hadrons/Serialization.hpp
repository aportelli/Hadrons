/*
 * Serialization.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
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

#ifndef Hadrons_Serialization_hpp_
#define Hadrons_Serialization_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

// ############# //
// FILE CONTENTS //
// ############# //
// This file consists of:
// ## Main Objects ##
// - GenericSerializable
//   This is a type-erased holder for Grid::Serializable subclasses and
//   collections thereof
//   Under a specific template specialisation for the reader and writer objects, 
//   GenericSerializables can be used to create collections of disparate
//   Grid::Serializables, e.g. std::vector<GenericSerializable<...>>.
//
// - SerializableGroup
//   This is a class that holds an std::vector of named GenericSerializables.
//
// ## General-purpose Typedefs ##
// - HadronsSerializable
//   A specialisation of GenericSerialzable with the default Hadrons reader
//   and writer classes.
//   NOTE: Most developers will want to use this typedef for type-erasing
//   Grid::Serializables.
//
// - HadronsSerializableGroup
//   A specialisation a SerializableGroup with the default Hadrons reader
//   and writer classes.
//   NOTE: Most developers will want to use this typedef for making collections
//   of Grid::Serializables.


// ############ //
// MAIN OBJECTS //
// ############ //

// type-erased holder for any Grid::Serializable descendant object
template<typename WriterType, typename ReaderType> // 
class GenericSerializable : Serializable
{
public:
    // ############ //
    // CONSTRUCTORS //
    // ############ //
    // Constructors
    GenericSerializable()             {}
    GenericSerializable(int i)        {} // Placeholder constructor that works with variadic macros
    template<typename T>
    GenericSerializable(const T& obj) { this->hold<T>(obj); }

    // Copy assignment constructors
    template<typename T>
    typename std::enable_if<!std::is_base_of<GenericSerializable, T>::value, GenericSerializable&>::type
    operator=(const T& obj)
    {
        this->hold<T>(obj);
        return *this;
    }

    // Get back object
    template<typename T>
    const T &get(void)
    {
        return dynamic_cast<const Model<T> *>(this->object.get())->serializable;
    }

    // ############ //
    // TYPE ERASURE //
    // ############ //
    template<typename T, typename... Args>
    T& hold(Args&&... args)
    {
        auto model_ptr = std::make_shared<Model<T>>(T{std::forward<Args>(args)...});
        this->object = model_ptr;
        return model_ptr->serializable;
    }

    struct Concept
    {
        virtual ~Concept() {}
        virtual void write(WriterType& wr, const std::string &s) const = 0;
    };

    template<typename U>
    struct Model final : Concept
    {
        U serializable;
        
        Model()                             {}
        Model(const U& u) : serializable(u) {}

        void write(WriterType& writer, const std::string &s) const override
        {
            writer.write(s, serializable);
        }
    };

    // ######### //
    // INTERFACE //
    // ######### //
    static inline void write(WriterType& wr, const std::string &s, const GenericSerializable& output)
    {
        output.object->write(wr, s);
    }

private:
    // ############### //
    // PRIVATE MEMBERS //
    // ############### //
    std::shared_ptr<const Concept> object=nullptr;
};

// a directly writeable holder for std::vectors of GenericSerializables that 
// ensures the group name is pushed onto the write stack when written
template<typename WriterType, typename ReaderType>
class SerializableGroup : Serializable
{
public:
    // ####### //
    // HELPERS //
    // ####### //
    typedef GenericSerializable<WriterType, ReaderType> Element_t;

    struct GroupElement
    {
        std::string name;
        Element_t serializable;
    };

    // ############# //
    // PULIC MEMBERS //
    // ############# //
    std::vector<GroupElement> elements;

    // ############ //
    // CONSTRUCTORS //
    // ############ //
    SerializableGroup()                         {}
    SerializableGroup(unsigned long elem_count) { this->elements.reserve(elem_count); }

    // ######### //
    // INTERFACE //
    // ######### //
    void append(const std::string& name, const Element_t& serializable)
    {
        this->elements.emplace_back(GroupElement{name, serializable});
    }

    template<typename T>
    void append(const std::string& name, const T& grid_serializable)
    {
        this->createElement<T>(name) = grid_serializable;
    }

    auto& createElement(const std::string& name)
    {
        this->elements.emplace_back(GroupElement{name, Element_t{}});
        auto& elem = this->elements.back();
        return elem.serializable;
    }

    template<typename T>
    T& createElement(const std::string& name)
    {
        return this->createElement(name).template hold<T>();
    }

    SerializableGroup& createSubGroup(const std::string& name)
    {
        return this->createElement<SerializableGroup>(name);
    }

    static inline void write(WriterType& wr, const std::string &s, const SerializableGroup& output)
    {
        if (!s.empty())
            wr.push(s);
        
        for (const auto& elem : output.elements)
        {
            auto& name = elem.name;
            auto& serializable = elem.serializable;
            serializable.write(wr, name, serializable);
        }   
        
        if (!s.empty())
            wr.pop();
    }
};

// ######################## //
// GENERAL-PURPOSE TYPEDEFS //
// ######################## //
typedef GenericSerializable<Grid::Writer<ResultWriter>, Grid::Writer<ResultReader>> HadronsSerializable;
typedef SerializableGroup  <Grid::Writer<ResultWriter>, Grid::Writer<ResultReader>> HadronsSerializableGroup;

END_HADRONS_NAMESPACE

#endif // Hadrons_Serialization_hpp_
