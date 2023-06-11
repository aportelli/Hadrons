/*
 * DistillationVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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
#ifndef Distillation_Vectors_hpp_
#define Distillation_Vectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  Methods for distillation vectors I/O                  *
 ******************************************************************************/
class DistillationVectorsIo
{
public:
    struct Record: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(Record,
                                        int, nNoise,
                                        int, nDL,
                                        int, nDS, 
                                        int, nDT, 
                                        std::vector<int>, timeSources,
                                        std::string, vecType,
                                        unsigned int, index);
        Record(void): index(0) {}
    };
public:
    template <typename Field>
    static void write(const std::string fileStem, 
                         std::vector<Field> &vec, 
                         const std::string vecType, 
                         const int nNoise, 
                         const int nDL,
                         const int nDS, 
                         const int nDT, 
                         std::vector<int> timeSources,
                         const bool multiFile, 
                         const int trajectory = -1);
    template <typename Field>
    static void read(std::vector<Field> &vec, 
                         const std::string fileStem,
                         const int nNoise, 
                         const int nDL,
                         const int nDS, 
                         const int nDT, 
                         const bool multiFile, 
                         const int trajectory = -1);
    template <typename Field>
    static void writeComponent(const std::string fileStem, 
                         Field &vec, 
                         const std::string vecType, 
                         const int nNoise, 
                         const int nDL,
                         const int nDS, 
                         const int nDT, 
                         std::vector<int> timeSources,
                         const int componentIndex, 
                         const int trajectory = -1);
    template <typename Field>
    static void readComponent(Field &vec, 
                         const std::string fileStem,
                         const int nNoise, 
                         const int nDL,
                         const int nDS, 
                         const int nDT, 
                         const int componentIndex, 
                         const int trajectory = -1);
private:
    static inline std::string vecFilename(const std::string stem, 
                                          const int traj, 
                                          const bool multiFile)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (multiFile)
        {
            return stem + t;
        }
        else
        {
            return stem + t + ".bin";
        }
    }
};


/******************************************************************************
 *               distillation vectors I/O template implementation               *
 ******************************************************************************/
template <typename Field>
void DistillationVectorsIo::write(const std::string fileStem, 
                                    std::vector<Field> &vec, 
                                    const std::string vecType, 
                                    const int nNoise, 
                                    const int nDL,
                                    const int nDS, 
                                    const int nDT, 
                                    std::vector<int> timeSources,
                                    const bool multiFile, 
                                    const int trajectory)
{
    Record       record;
    GridBase     *grid = vec[0].Grid();
    ScidacWriter binWriter(grid->IsBoss());
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    record.vecType = vecType;
    record.nNoise = nNoise;
    record.nDL = nDL;
    record.nDS = nDS;
    record.nDT = nDT;
    record.timeSources = timeSources;
    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Writing vector " << i << std::endl;
            makeFileDir(fullFilename, grid);
            binWriter.open(fullFilename);
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
            binWriter.close();
        }
    }
    else
    {
        makeFileDir(filename, grid);
        binWriter.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Writing vector " << i << std::endl;
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
        }
        binWriter.close();
    }
}

template <typename Field>
void DistillationVectorsIo::read(std::vector<Field> &vec, 
                                    const std::string fileStem, 
                                    const int nNoise, 
                                    const int nDL,
                                    const int nDS, 
                                    const int nDT, 
                                    const bool multiFile, 
                                    const int trajectory)
{
    Record       record;
    ScidacReader binReader;
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.open(fullFilename);
            binReader.readScidacFieldRecord(vec[i], record);
            binReader.close();
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
            if (record.nNoise != nNoise || record.nDL != nDL || record.nDS != nDS || record.nDT != nDT )
            {
                HADRONS_ERROR(Io, "dilution parameter mismatch");
            }
        }
    }
    else
    {
        binReader.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.readScidacFieldRecord(vec[i], record);
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
        binReader.close();
    }
}


/******************************************************************************
 *               distillation vectors I/O template implementation             *
 *               version for single component                                 *
 ******************************************************************************/
template <typename Field>
void DistillationVectorsIo::writeComponent(const std::string fileStem, 
                                    Field &vec, 
                                    const std::string vecType, 
                                    const int nNoise, 
                                    const int nDL,
                                    const int nDS, 
                                    const int nDT, 
                                    std::vector<int> timeSources,
                                    const int componentIndex, 
                                    const int trajectory)
{
    Record       record;
    GridBase     *grid = vec.Grid();
    ScidacWriter binWriter(grid->IsBoss());
    std::string  filename = vecFilename(fileStem, trajectory, 1);

    record.vecType = vecType;
    record.nNoise = nNoise;
    record.nDL = nDL;
    record.nDS = nDS;
    record.nDT = nDT;
    record.timeSources = timeSources;
    std::string fullFilename;

    fullFilename = filename + "/elem" + std::to_string(componentIndex) + ".bin";

    LOG(Message) << "Writing vector " << componentIndex << std::endl;
    makeFileDir(fullFilename, grid);
    binWriter.open(fullFilename);
    record.index = componentIndex;
    binWriter.writeScidacFieldRecord(vec, record);
    binWriter.close();
}

template <typename Field>
void DistillationVectorsIo::readComponent(Field &vec, 
                                    const std::string fileStem, 
                                    const int nNoise, 
                                    const int nDL,
                                    const int nDS, 
                                    const int nDT, 
                                    const int componentIndex, 
                                    const int trajectory)
{
    Record       record;
    ScidacReader binReader;
    std::string  filename = vecFilename(fileStem, trajectory, 1);

    std::string fullFilename;

    fullFilename = filename + "/elem" + std::to_string(componentIndex) + ".bin";

    LOG(Message) << "Reading vector " << componentIndex << std::endl;
    binReader.open(fullFilename);
    binReader.readScidacFieldRecord(vec, record);
    binReader.close();
    if (record.index != componentIndex)
    {
        HADRONS_ERROR(Io, "vector index mismatch");
    }
    if (record.nNoise != nNoise || record.nDL != nDL || record.nDS != nDS || record.nDT != nDT )
    {
        HADRONS_ERROR(Io, "dilution parameter mismatch");
    }
}
END_HADRONS_NAMESPACE

#endif // Distillation_Vectors_hpp_
