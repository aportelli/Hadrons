/*
 * SaveNersc.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_SaveNersc_hpp_
#define Hadrons_MIO_SaveNersc_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Save a NERSC configuration                           *
 ******************************************************************************/

// gauge          Name of the gauge field object to write
// file           Name of the file to write the gauge field to
// gaugeMetadata  (optional) Name of the file the gauge was loaded from
//                ... so metadata can be copied
// comment        (optional) value to write in COMMENT field of gauge-field metadata
//                e.g. "gauge fixed: alpha=0.05, maxiter=1000000, Omega_tol=1e-12,
//                      Phi_tol=1e-12, gaugeFix=coulomb, Fourier=true"

BEGIN_MODULE_NAMESPACE(MIO)

class SaveNerscPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveNerscPar,
                                    std::string, gauge,
                                    std::string, file,
                                    std::string, gaugeMetadata,
                                    std::string, comment);
};

template <typename GImpl>
class TSaveNersc: public Module<SaveNerscPar>
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TSaveNersc(const std::string name);
    // destructor
    virtual ~TSaveNersc(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveNersc,  TSaveNersc<GIMPL>,  MIO);

/******************************************************************************
*                       TSaveNersc implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TSaveNersc<GImpl>::TSaveNersc(const std::string name)
: Module<SaveNerscPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TSaveNersc<GImpl>::getInput(void)
{
  return { par().gauge };
}

template <typename GImpl>
std::vector<std::string> TSaveNersc<GImpl>::getOutput(void)
{
    return {};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl>
void TSaveNersc<GImpl>::setup(void)
{
}

// Class to write metadata /////////////////////////////////////////////////////
template<typename GImpl> class CachedGaugeFixedStatistics
{
public:
    GAUGE_TYPE_ALIASES(GImpl,);
    FieldMetaData fmd;

    CachedGaugeFixedStatistics(std::string metadataFile, GridBase *grid)
    {
        LOG(Message) << "Copying NERSC metadata from file '" << metadataFile
                     << "'" << std::endl;
        NerscIO::readHeader(metadataFile, grid, fmd);
    }

    void operator()(GaugeField & data,FieldMetaData &header)
    {
        header.link_trace=WilsonLoops<GImpl>::linkTrace(data);
        header.plaquette =WilsonLoops<GImpl>::avgPlaquette(data);
        header.hdr_version = fmd.hdr_version;
        header.storage_format = fmd.storage_format;
        header.ensemble_id = fmd.ensemble_id;
        header.ensemble_label = fmd.ensemble_label;
        header.sequence_number = fmd.sequence_number;
    }
};

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl>
void TSaveNersc<GImpl>::execute(void)
{
    std::map<std::string,std::string> extraMetadata;
    if( !par().comment.empty() )
    {
        extraMetadata["COMMENT"] = par().comment;
    }

    std::string   traj = "." + std::to_string(vm().getTrajectory());
    std::string fileName = par().file + traj;
    LOG(Message) << "Saving NERSC configuration to file '" << fileName
                 << "'" << std::endl;

    auto &U = envGet(GaugeField, par().gauge);
    const bool bCopyMetadata = !par().gaugeMetadata.empty();
    if( bCopyMetadata )
    {
        using CachedStats = CachedGaugeFixedStatistics<GImpl>;
        std::string metadataFile = par().gaugeMetadata + traj;
        CachedStats Cache(metadataFile, U.Grid());
        NerscIO::writeConfiguration<CachedGaugeFixedStatistics<GImpl>>
                                   (U, fileName, 0, 0, Cache, &extraMetadata);
    }
    else
    {
        NerscIO::writeConfiguration<GaugeStatistics<GImpl>>(U, fileName, 0, 0, &extraMetadata);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_SaveNersc_hpp_
