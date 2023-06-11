/*
 * NamedTensor.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 *  Author: Felix Erben <ferben@ed.ac.uk>
 *  Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: ferben <ferben@debian.felix.com>
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

#ifndef Hadrons_NamedTensor_hpp_
#define Hadrons_NamedTensor_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 NamedTensor contains:
 1) Name of the tensor. By default, this is the tag name used for save / load
 2) Eigen::Tensor of type Scalar_ and rank NumIndices_ (row-major order)
 3) Name for each index
 They can be persisted to / restored from disk. During restore, these validations are performed:
   1) Tensor dimensionality must match
   2) IndexNames are validated against current values
   3) If the tensor has non-zero size, the tensor being loaded must have same extent in each dimension
 ******************************************************************************/

class NamedTensorDefaultMetadata : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NamedTensorDefaultMetadata, std::string, Version);
};

template<typename Scalar_, int NumIndices_, typename MetaData_ = NamedTensorDefaultMetadata>
class NamedTensor : Serializable
{
public:
    using Scalar = Scalar_;
    static constexpr int NumIndices = NumIndices_;
    using ET = Eigen::Tensor<Scalar_, NumIndices_, Eigen::RowMajor>;
    using Index = typename ET::Index;
    using Traits = Grid::EigenIO::Traits<ET>;
    protected:
      Grid::Vector<typename Traits::scalar_type> deviceBuf;
    public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(NamedTensor,
                                    Eigen::TensorMap<ET>,     tensor,
                                    std::vector<std::string>, IndexNames,
                                    MetaData_,                MetaData );

    // Name of the object and Index names as set in the constructor
    const std::string                          &Name_;
    const std::array<std::string, NumIndices_> &DefaultIndexNames_;

    // Return the product of a variable number of dimensions
    Eigen::Index NamedTensorSize() { return 1; }
    template<typename... IndexTypes>
    inline Eigen::Index NamedTensorSize(Eigen::Index firstDimension, IndexTypes... otherDimensions)
    {
        if( sizeof...(otherDimensions) )
            firstDimension *= NamedTensorSize( otherDimensions... );
        return firstDimension;
    }

    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    NamedTensor(const std::string &Name,
                                                      const std::array<std::string, NumIndices_> &indexNames,
                                                      Eigen::Index firstDimension, IndexTypes... otherDimensions)
    : deviceBuf( Traits::count * NamedTensorSize(firstDimension, otherDimensions...) ),
      tensor(reinterpret_cast<Scalar *>(&deviceBuf[0]), firstDimension, otherDimensions...),
      IndexNames{indexNames.begin(), indexNames.end()}, Name_{Name}, DefaultIndexNames_{indexNames}
    {
        if(sizeof...(otherDimensions) + 1 != NumIndices_)
        {
            HADRONS_ERROR(Argument, "NamedTensor: dimensions != tensor rank");
        }
    }

    // Do my index names match the default for my type?
    template<typename array_or_vector_of_string>
    bool ValidateIndexNames( const array_or_vector_of_string &CheckNames ) const
    {
        return IndexNames.size() == CheckNames.size() && std::equal( IndexNames.begin(), IndexNames.end(), CheckNames.begin(),
            [](const std::string &s1, const std::string &s2)
            {
                 return s1.size() == s2.size() && std::equal( s1.begin(), s1.end(), s2.begin(),
                     [](const char & c1, const char & c2)
                     { return c1 == c2 || std::toupper(c1) == std::toupper(c2); }); // case insensitive
            });
    }
    bool ValidateIndexNames() const { return ValidateIndexNames(DefaultIndexNames_); }

    void write(const std::string &FileName, const std::string &Tag) const
    {
        #ifdef HAVE_HDF5
        std::string FileName_{FileName};
        FileName_.append( ".h5" );
        LOG(Message) << "Writing " << Name_ << " to file " << FileName_ << " tag " << Tag << std::endl;
        using ScalarType = typename Traits::scalar_type;        
        std::vector<hsize_t> dims, 
                         gridDims;

        constexpr unsigned int ContainerRank{Traits::Rank};     
        // These are the tensor dimensions
        for (int i = 0; i < NumIndices_; i++)
        {
            dims.push_back(tensor.dimension(i));
        }
        // These are the dimensions of the Grid datatype - for vectors and matrices, we do not write dimensions of size 1
        for (int i = 0; i < ContainerRank; i++)
        {
            if(Traits::Dimension(i) > 1)
            {
                dims.push_back(Traits::Dimension(i));
            }
            gridDims.push_back(Traits::Dimension(i));
        }  
        //Convention for scalar containers - without this the writing of gridDims fails 
        if(gridDims.empty())
        {
            gridDims.push_back(1);
        }        
    
        Hdf5Writer writer( FileName_ );
        Grid::write (writer, "MetaData", MetaData);
        Grid::write (writer, "IndexNames", IndexNames);
        Grid::write (writer, "GridDimensions", gridDims);
        Grid::write (writer, "TensorDimensions", dims);
        H5NS::DataSet dataset;
        H5NS::DataSpace      dataspace(dims.size(), dims.data());
        H5NS::DSetCreatPropList     plist;
        
        plist.setFletcher32();
        plist.setChunk(dims.size(), dims.data());
        H5NS::Group &group = writer.getGroup();
        dataset     = group.createDataSet(Tag,Hdf5Type<ScalarType>::type(), dataspace, plist);

        dataset.write(tensor.data(),Hdf5Type<ScalarType>::type(), dataspace);
        #else
        HADRONS_ERROR(Implementation, "NamedTensor I/O needs HDF5 library");
        #endif
    }
    void write(const std::string &FileName) const { return write(FileName, Name_); }

    // Read tensor.
    // Validate:
    //  1) index names (if requested)
    //  2) index dimensions (if they are non-zero when called)
    template<typename Reader> void read(Reader &reader, bool bValidate, const std::string &Tag)
    {
        #ifdef HAVE_HDF5
        // Grab index names and dimensions
        std::vector<std::string> OldIndexNames{std::move(IndexNames)};
        const typename ET::Dimensions OldDimensions{tensor.dimensions()};

        using ScalarType = typename Traits::scalar_type;        
        std::vector<hsize_t> dims,
                         gridDims;
    
        Grid::read (reader, "GridDimensions", gridDims);
        Grid::read (reader, "TensorDimensions", dims);
        Grid::read (reader, "MetaData", MetaData);
        Grid::read (reader, "IndexNames", IndexNames);
    
        H5NS::DataSet dataset;
        H5NS::Group &group = reader.getGroup();
        dataset=group.openDataSet(Tag);
        dataset.read(tensor.data(),Hdf5Type<ScalarType>::type());

        //validate dimesnions and labels
        const typename ET::Dimensions & NewDimensions{tensor.dimensions()};
        for (int i = 0; i < NumIndices_; i++)
            if(OldDimensions[i] && OldDimensions[i] != NewDimensions[i])
            {
                HADRONS_ERROR(Size,"NamedTensor::read dimension size");
            }
        if (bValidate && !ValidateIndexNames(OldIndexNames))
        {
            HADRONS_ERROR(Definition,"NamedTensor::read dimension name");
        }
        #else
        HADRONS_ERROR(Implementation, "NamedTensor I/O needs HDF5 library");
        #endif
    }
    template<typename Reader> void read(Reader &r, bool bValidate = true) { read(r, bValidate, Name_); }

    inline void read (const std::string &FileName, bool bValidate, const std::string &Tag)
    {
        #ifdef HAVE_HDF5
        std::string FileName_{FileName};
        FileName_.append( ".h5" );
        LOG(Message) << "reading " << FileName_ << std::endl;
        Hdf5Reader r( FileName_ );
        read(r, bValidate, Tag);
        #else
        HADRONS_ERROR(Implementation, "NamedTensor I/O needs HDF5 library");
        #endif
    }
    inline void read (const std::string &FileName, bool bValidate= true) { return read(FileName, bValidate, Name_); }
};

/******************************************************************************
 Common elements for distillation
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MDistil)

class PerambMetadata : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PerambMetadata, 
                                    std::string, Version,
                                    std::vector<int>, timeSources );
};

class PerambTensor : public NamedTensor<SpinVector, 6, PerambMetadata>
{
    public:
    static const std::string                Name__;
    static const std::array<std::string, 6> DefaultIndexNames__;
    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    PerambTensor(Eigen::Index nT, Eigen::Index nVec, Eigen::Index nDl, Eigen::Index nNoise, Eigen::Index nTinv, Eigen::Index nDs)
    : NamedTensor{Name__, DefaultIndexNames__, nT, nVec, nDl, nNoise, nTinv, nDs} {}
};

// Separate class for multiFile
class PerambIndexMetadata : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PerambIndexMetadata, 
                                    std::string, Version,
                                    int, timeDilutionIndex,
                                    std::vector<std::string>, noiseHashes );
};

class PerambIndexTensor : public NamedTensor<SpinVector, 5, PerambIndexMetadata>
{
    public:
    static const std::string                Name__;
    static const std::array<std::string, 5> DefaultIndexNames__;
    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    PerambIndexTensor(Eigen::Index nT, Eigen::Index nVec, Eigen::Index nDl, Eigen::Index nNoise, Eigen::Index nDs)
    : NamedTensor{Name__, DefaultIndexNames__, nT, nVec, nDl, nNoise, nDs} {}
};

class TimesliceEvals : public NamedTensor<RealD, 2>
{
    public:
    static const std::string                Name__;
    static const std::array<std::string, 2> DefaultIndexNames__;
    // Construct a named tensor explicitly specifying size of each dimension
    template<typename... IndexTypes>
    TimesliceEvals(Eigen::Index nT, Eigen::Index nVec)
    : NamedTensor{Name__, DefaultIndexNames__, nT, nVec} {}
};

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_NamedTensor_hpp_
