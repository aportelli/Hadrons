/*
 * FieldIo.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_FieldIo_hpp_
#define Hadrons_FieldIo_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE


/******************************************************************************
 *                      Multi-precision field writer                          *
 ******************************************************************************/
template <typename T, typename TIo = T>
class FieldWriter
{
public:
    // constructors
    FieldWriter(void) = default;
    FieldWriter(GridBase *grid, GridBase *gridIo = nullptr);
    FieldWriter(const std::string filename, GridBase *grid, GridBase *gridIo = nullptr);
    // destructor
    virtual ~FieldWriter(void);
    // set grids
    void setGrids(GridBase *grid, GridBase *gridIo = nullptr);
    // open/close file
    void open(const std::string filename);
    void open(const std::string filename, GridBase *grid, GridBase *gridIo = nullptr);
    void close(void);
    // write data
    void writeHeader(XmlWriter &xmlWriter, const std::string name);
    template <typename Metadata>
    void writeField(T &field, const Metadata &md);
private:
    std::string                   filename_;
    bool                          isOpened_{false};
    GridBase                      *grid_, *gridIo_{nullptr};
    std::unique_ptr<T>            testBuf_{nullptr};
    std::unique_ptr<TIo>          ioBuf_{nullptr};
    std::unique_ptr<ScidacWriter> binWriter_{nullptr};
};

/******************************************************************************
 *                      Multi-precision field reader                          *
 ******************************************************************************/
template <typename T, typename TIo = T>
class FieldReader
{
public:
    // constructors
    FieldReader(void) = default;
    FieldReader(GridBase *grid, GridBase *gridIo = nullptr);
    FieldReader(const std::string filename, GridBase *grid, GridBase *gridIo = nullptr);
    // destructor
    virtual ~FieldReader(void);
    // set grids
    void setGrids(GridBase *grid, GridBase *gridIo = nullptr);
    // open/close file
    void open(const std::string filename);
    void open(const std::string filename, GridBase *grid, GridBase *gridIo = nullptr);
    void close(void);
    // read data
    void readHeader(std::string &xmlString);
    template <typename Metadata>
    void readField(T &field, Metadata &md);
private:
    std::string                   filename_;
    bool                          isOpened_{false};
    GridBase                      *grid_, *gridIo_{nullptr};
    std::unique_ptr<TIo>          ioBuf_{nullptr};
    std::unique_ptr<ScidacReader> binReader_{nullptr};
};

/******************************************************************************
 *                      FieldWriter template implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T, typename TIo>
FieldWriter<T, TIo>::FieldWriter(GridBase *grid, GridBase *gridIo)
{
    setGrids(grid, gridIo);
}

template <typename T, typename TIo>
FieldWriter<T, TIo>::FieldWriter(const std::string filename, GridBase *grid, GridBase *gridIo)
{
    open(filename, grid, gridIo);
}

// destructor //////////////////////////////////////////////////////////////////
template <typename T, typename TIo>
FieldWriter<T, TIo>::~FieldWriter(void)
{
    close();
}

// set grids ///////////////////////////////////////////////////////////////////
template <typename T, typename TIo>
void FieldWriter<T, TIo>::setGrids(GridBase *grid, GridBase *gridIo)
{
    grid_   = grid;
    gridIo_ = gridIo;
    if (typeHash<T>() != typeHash<TIo>())
    {
        if (gridIo_ == nullptr)
        {
            HADRONS_ERROR(Definition, 
                          "I/O type different from field type but null I/O grid passed");
        }
        ioBuf_.reset(new TIo(gridIo_));
        testBuf_.reset(new T(grid_));
    }
}

// open/close file /////////////////////////////////////////////////////////////
template <typename T, typename TIo>
void FieldWriter<T, TIo>::open(const std::string filename)
{
    if (grid_ == nullptr)
    {
        HADRONS_ERROR(Definition, "no grid has been set");
    }
    if (isOpened_)
    {
        if (filename != filename_)
        {
            close();
        }
        else
        {
            return;
        }
    }
    filename_ = filename;
    makeFileDir(filename_, grid_);
    binWriter_.reset(new ScidacWriter(grid_->IsBoss()));
    binWriter_->open(filename_);
    isOpened_ = true;
}

template <typename T, typename TIo>
void FieldWriter<T, TIo>::open(const std::string filename, GridBase *grid, GridBase *gridIo)
{
    setGrids(grid, gridIo);
    open(filename);
}

template <typename T, typename TIo>
void FieldWriter<T, TIo>::close(void)
{
    if (isOpened_)
    {
        filename_ = "";
        binWriter_->close();
        binWriter_.reset(nullptr);
        isOpened_ = false;
    }
}

// write data //////////////////////////////////////////////////////////////////
template <typename T, typename TIo>
void FieldWriter<T, TIo>::writeHeader(XmlWriter &xmlWriter, const std::string name)
{
    if (isOpened_)
    {
        binWriter_->writeLimeObject(1, 1, xmlWriter, name, SCIDAC_FILE_XML);
    }
    else
    {
        HADRONS_ERROR(Io, "No file opened");
    }
}

template <typename T, typename TIo>
template <typename Metadata>
void FieldWriter<T, TIo>::writeField(T &field, const Metadata &md)
{
    if (isOpened_)
    {
        if (typeHash<T>() == typeHash<TIo>())
        {
            binWriter_->writeScidacFieldRecord(field, md, DEFAULT_ASCII_PREC);
        }
        else
        {
            precisionChange(*ioBuf_, field);
            precisionChange(*testBuf_, *ioBuf_);
            *testBuf_ -= field;
            LOG(Message) << "Precision diff norm^2 " << norm2(*testBuf_) << std::endl;
            binWriter_->writeScidacFieldRecord(*ioBuf_, md, DEFAULT_ASCII_PREC);
        }
    }
    else
    {
        HADRONS_ERROR(Io, "No file opened");
    }
}

/******************************************************************************
 *                      FieldReader template implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template<typename T, typename TIo>
FieldReader<T, TIo>::FieldReader(GridBase *grid, GridBase *gridIo)
{
    setGrids(grid, gridIo);
}

template<typename T, typename TIo>
FieldReader<T, TIo>::FieldReader(const std::string filename, 
                                 GridBase *grid, GridBase *gridIo)
{
    open(filename, grid, gridIo);
}

// destructor //////////////////////////////////////////////////////////////////
template<typename T, typename TIo>
FieldReader<T, TIo>::~FieldReader(void)
{
    close();
}

// set grids ///////////////////////////////////////////////////////////////////
template <typename T, typename TIo>
void FieldReader<T, TIo>::setGrids(GridBase *grid, GridBase *gridIo)
{
    grid_   = grid;
    gridIo_ = gridIo;
    if (typeHash<T>() != typeHash<TIo>())
    {
        if (gridIo_ == nullptr)
        {
            HADRONS_ERROR(Definition, 
                          "I/O type different from field type but null I/O grid passed");
        }
        ioBuf_.reset(new TIo(gridIo_));
    }
}

// open/close file /////////////////////////////////////////////////////////////
template<typename T, typename TIo>
void FieldReader<T, TIo>::open(const std::string filename)
{
    if (isOpened_)
    {
        if (filename != filename_)
        {
            close();
        }
        else
        {
            return;
        }
    }
    filename_ = filename;
    binReader_.reset(new ScidacReader);
    binReader_->open(filename_);
    isOpened_ = true;
}

template<typename T, typename TIo>
void FieldReader<T, TIo>::open(const std::string filename, 
                               GridBase *grid, GridBase *gridIo)
{
    setGrids(grid, gridIo);
    open(filename);
}

template<typename T, typename TIo>
void FieldReader<T, TIo>::close(void)
{
    if (isOpened_)
    {
        filename_ = "";
        binReader_->close();
        binReader_.reset(nullptr);
        isOpened_ = false;
    }
}

// read data ///////////////////////////////////////////////////////////////////
template<typename T, typename TIo>
void FieldReader<T, TIo>::readHeader(std::string &xmlString)
{
    if (isOpened_)
    {
        binReader_->readLimeObject(xmlString, SCIDAC_FILE_XML);
    }
    else
    {
        HADRONS_ERROR(Io, "No file opened");
    }
}

template<typename T, typename TIo>
template <typename Metadata>
void FieldReader<T, TIo>::readField(T &field, Metadata &md)
{
    if (isOpened_)
    {
        if (typeHash<T>() == typeHash<TIo>())
        {
            binReader_->readScidacFieldRecord(field, md);
        }
        else
        {
            binReader_->readScidacFieldRecord(*ioBuf_, md);
            precisionChange(field, *ioBuf_);
        }
    }
    else
    {
        HADRONS_ERROR(Io, "No file opened");
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_FieldIo_hpp_
