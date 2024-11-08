// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCOLLISIONSDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCOLLISIONSDERIVED_H_

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"

namespace o2::aod
{
namespace femtocollisions
{
DECLARE_SOA_COLUMN(Mult, mult, float);             //! Multiplicity estimator set by producer
DECLARE_SOA_COLUMN(Cent, cent, float);             //! Centrality estimator (or multiplicity percentile estimator) set by producer
DECLARE_SOA_COLUMN(MagField, magField, float);     //! Magnetic field of the event
DECLARE_SOA_COLUMN(Sphericity, sphericity, float); //! Sphericity of the event
} // namespace femtocollisions

// table for basic collision information
DECLARE_SOA_TABLE_VERSIONED(FCols_001, "AOD", "FCOLS", 1,
                            o2::soa::Index<>,           //! Index
                            o2::aod::collision::PosZ,   //! z coordinate of vertex
                            femtocollisions::Mult,      //! multiplicity
                            femtocollisions::Cent,      //! centrality
                            femtocollisions::MagField); //! magnetic field
using FCols = FCols_001;

// table for for primary vertex location
DECLARE_SOA_TABLE_VERSIONED(FColPos_001, "AOD", "FCOLPOS", 1,
                            collision::PosX,  //! x coordinate of vertex
                            collision::PosY); //! y coordindate of vertex
using FColPos = FColPos_001;

// table for different multiplicity estimators
DECLARE_SOA_TABLE_VERSIONED(FColMults_001, "AOD", "FCOLMULTS", 1,
                            evsel::Sel8, evsel::Selection,                  //! event selection: sel8
                            mult::MultFT0A, mult::MultFT0C, mult::MultFV0A, //! FIT detectors
                            mult::MultNTracksPVeta1,                        //! track multiplicities with eta cut for INEL>0
                            mult::MultPVTotalContributors,                  //! number of PV contribs total
                            mult::MultNTracksGlobal,                        //! global track multiplicities
                            mult::MultNTracksITSTPC,                        //! track multiplicities, PV contribs, no eta cut
                            mult::MultAllTracksTPCOnly,                     //! TPConly track multiplicities, all, no eta cut
                            mult::MultAllTracksITSTPC,                      //! ITSTPC track multiplicities, all, no eta cut
                            mult::MultZNA, mult::MultZNC, mult::MultZEM1,   //! ZDC signals
                            mult::MultZEM2, mult::MultZPA, mult::MultZPC,
                            evsel::NumTracksInTimeRange); //! add occupancy as extra
using FColMults = FColMults_001;

// table for different centrality (multiplicity percentile) estimators
DECLARE_SOA_TABLE_VERSIONED(FColCents_001, "AOD", "FCOLCENTS", 1,
                            cent::CentFT0M, //! centrality from FT0M
                            cent::CentFT0A, //! centrality from FT0A
                            cent::CentFT0C, //! centrality from FT0c
                            cent::CentFV0A);
using FColCents = FColCents_001;

// table for collision sphericity
DECLARE_SOA_TABLE_VERSIONED(FColSph_001, "AOD", "FCOLSPH", 1,
                            femtocollisions::Sphericity); //! collision sphericity
using FColSph = FColSph_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCOLLISIONSDERIVED_H_
