/*
 * FieldIoSharded.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <Ryan.Hill@ed.ac.uk>
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
#ifndef Hadrons_FieldIoSharded_hpp_
#define Hadrons_FieldIoSharded_hpp_

#include <Hadrons/Global.hpp>


BEGIN_HADRONS_NAMESPACE


static constexpr inline uint32_t FIELDIOSHARDED_MAGIC_VALUE = 0x53485244; // 'SHRD'

// alignas() in case we try O_DIRECT in the future (0x4000 & different flags for Darwin)
struct alignas(0x1000) FieldIoShardedFileHeader
{
    uint32_t magic;
    uint32_t payloadCrc;
    uint64_t payloadBytes;
    int32_t  localStarts[4];
    int32_t  localDimensions[4];
    int32_t  globalDimensions[4];
};

/**********************************************************************
 *                      Sharded Field Writer                          *
 **********************************************************************/
template <typename T>
void writeShardedFile(const std::string& filepath, const T& field)
{
    // Fill in data header
    GridBase *grid = field.Grid();
    FieldIoShardedFileHeader hdr;
    std::memset(&hdr, 0, sizeof(hdr));
    hdr.magic = FIELDIOSHARDED_MAGIC_VALUE;
    size_t size, sizec;
    uint32_t crc;
    GridStopWatch ioWatch, crcWatch;

    size = field.Grid()->oSites() * sizeof(typename T::vector_object);
    sizec = size / sizeof(char);

    autoView(field_v, field, CpuRead);
    crcWatch.Start();
    hdr.payloadCrc = GridChecksum::crc32(field_v.cpu_ptr, size);
    crcWatch.Stop();
    hdr.payloadBytes = sizec;
    for (int i=0; i < 4; ++i)
    {
        hdr.localStarts[i]      = field.Grid()->LocalStarts()[i];
        hdr.localDimensions[i]  = field.Grid()->LocalDimensions()[i];
        hdr.globalDimensions[i] = field.Grid()->GlobalDimensions()[i];
    }
    
    // Data write
    ioWatch.Start();
    std::ofstream file(filepath, std::ios::out | std::ios::binary);
    if (!file.is_open())
    {
        HADRONS_ERROR(Io, "Failed to open file " + filepath);
    }
    if (!file.write(reinterpret_cast<char *>(&hdr), sizeof(hdr)))
    {
        HADRONS_ERROR(Io, "Failed to write header ");
    }
    if(!file.write(reinterpret_cast<char *>(field_v.cpu_ptr), sizec))
    {
        HADRONS_ERROR(Io, "Failed to write payload ");
    }
    if (!file.flush())
    {
        HADRONS_ERROR(Io, "Failed to flush file");
    }
    assert(file.tellp() == (sizeof(hdr) + sizec));
    file.close();
    if (!file)
    {
        HADRONS_ERROR(Io, "Failed to close file");
    }
    
    ioWatch.Stop();
    
    // Reporting
    size *= field.Grid()->ProcessorCount();
    auto time = ioWatch.useconds();
    auto mbytesPerSecond = size / 1024. / 1024. / (ioWatch.useconds() / 1.e6);

    LOG(Message) << "Wrote " << size << " bytes in " << ioWatch.Elapsed()
                 << ", " << mbytesPerSecond << " MB/s" << std::endl;
    LOG(Message) << "Checksum overhead " << crcWatch.Elapsed() << std::endl;
}

/**********************************************************************
 *                      Sharded Field Reader                          *
 **********************************************************************/
template <typename T>
void readShardedFile(const std::string& filepath, T& field)
{
    // Data read
    GridBase *grid = field.Grid();
    FieldIoShardedFileHeader hdr;
    size_t size, sizec;
    uint32_t crcData;
    GridStopWatch ioWatch, crcWatch;

    size = field.Grid()->oSites() * sizeof(typename T::vector_object);
    sizec = size / sizeof(char);

    std::ifstream file(filepath, std::ios::in | std::ios::binary);
    if (!file.is_open())
    {
        HADRONS_ERROR(Io, "Failed to open file " + filepath);
    }

    // Read
    {
        autoView(field_v, field, CpuWrite);
        ioWatch.Start();
        if(!file.read(reinterpret_cast<char *>(&hdr), sizeof(hdr)))
        {
            HADRONS_ERROR(Io, "Failed to read header");
        }
        
        if (sizec != hdr.payloadBytes)
        {        
            std::stringstream err_str;
            err_str << "Unexpected payload size, "
                    << ": file header " << hdr.payloadBytes
                    << ", field size "     << sizec;
        }
        
        if(!file.read(reinterpret_cast<char *>(field_v.cpu_ptr), sizec))
        {
            HADRONS_ERROR(Io, "Failed to read payload");
        }
        
        assert(file.tellg() == sizeof(hdr) + sizec);
        ioWatch.Stop();
        file.close();
    }
    
    {
      autoView(field_v, field, CpuRead);
      crcWatch.Start();
      crcData = GridChecksum::crc32(field_v.cpu_ptr, size);
      crcWatch.Stop();
    }

    // Validation
    if (hdr.magic != FIELDIOSHARDED_MAGIC_VALUE) // SHRD
    {
        std::stringstream err_str;
        err_str << "Not a 'SHRD' file: magic value is 0x"
                << std::hex << hdr.magic << std::dec << " (expected 0x"
                << std::hex << FIELDIOSHARDED_MAGIC_VALUE << std::dec << ")";
        HADRONS_ERROR(Io, err_str.str());        
    }
    for (int i=0; i < 4; ++i)
    {
        if(hdr.localStarts[i]      != field.Grid()->LocalStarts()[i])
        {
            int32_t* a = &hdr.localStarts[0];
            std::stringstream err_str;
            err_str << "LocalStarts does not match runtime value: "
                    << "file header [" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << "], "
                    << "runtime " << field.Grid()->LocalStarts();
            HADRONS_ERROR(Io, err_str.str());
        }
        if(hdr.localDimensions[i]  != field.Grid()->LocalDimensions()[i])
        {
            int32_t* a = &hdr.localDimensions[0];
            std::stringstream err_str;
            err_str << "LocalDimensions does not match runtime value: "
                    << "file header [" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << "], "
                    << "runtime " << field.Grid()->LocalDimensions() << "";
            HADRONS_ERROR(Io, err_str.str());
        }
        if(hdr.globalDimensions[i] != field.Grid()->GlobalDimensions()[i])
        {
            int32_t* a = &hdr.globalDimensions[0];
            std::stringstream err_str;
            err_str << "GlobalDimensions does not match runtime value: "
                    << "file header [" << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << "], "
                    << "runtime " << field.Grid()->GlobalDimensions() << "";
            HADRONS_ERROR(Io, err_str.str());
        }
    }
    if(crcData != hdr.payloadCrc)
    {
        std::stringstream err_str;
        err_str << "File CRC does not match on rank " << grid->ThisRank()
                << ": file header 0x" << std::hex << hdr.payloadCrc << std::dec
                << ", payload 0x"     << std::hex << crcData        << std::dec;              
        HADRONS_ERROR(Io, err_str.str());
    }

    // Reporting
    size *= field.Grid()->ProcessorCount();
    auto time = ioWatch.useconds();
    auto mbytesPerSecond = size / 1024. / 1024. / (ioWatch.useconds() / 1.e6);

    LOG(Message) << "Read " << size << " bytes in " << ioWatch.Elapsed()
                 << ", " << mbytesPerSecond << " MB/s" << std::endl;
    LOG(Message) << "Checksum overhead " << crcWatch.Elapsed() << std::endl;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_FieldIoSharded_hpp_
