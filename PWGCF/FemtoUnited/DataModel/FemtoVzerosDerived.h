// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOVZEROSDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOVZEROSDERIVED_H_

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

namespace o2::aod
{
namespace femtovzeros
{

// columns for Vzero
DECLARE_SOA_COLUMN(VzeroMass, vzeroMass, float);                         //! Mass of Vzero
DECLARE_SOA_COLUMN(VzeroMask, vzeroMask, femtodatatypes::VzeroMaskType); //! Bitmask for Vzero selections

// columns for Vzero debug information
DECLARE_SOA_COLUMN(DauDCA, dauDCA, float);           //! Vzero daughter DCA at decay vertex
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Vzero transvers radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! x coordinate of Vzero decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! y coordinate of Vzero decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! z coordinate of Vzero decay vertex
DECLARE_SOA_COLUMN(KaonMass, kaonMass, float);       //! Vzero mass using Kaon hypothesis

// columns for Vzero daughter tracks
DECLARE_SOA_INDEX_COLUMN(Track, track);                                                //! Index in track table
DECLARE_SOA_COLUMN(DauTrackMask, dauTrackMask, femtodatatypes::VzeroDauTrackMaskType); //! Bitmask for Vzero daughter track selections
DECLARE_SOA_COLUMN(DauTPCMask, dauTPCMask, femtodatatypes::VzeroDauTPCMaskType);       //! Bitmask for Vzero daughter TPC PID selections

} // namespace femtovzeros

// table for basic vzero information
DECLARE_SOA_TABLE_VERSIONED(FVzeros_001, "AOD", "FVZEROS", 1,
                            o2::soa::Index<>,
                            femtobase::CollisionId,
                            femtobase::Pt,
                            femtobase::Eta,
                            femtobase::Phi,
                            femtovzeros::VzeroMass,
                            femtobase::Theta<femtobase::Eta>,
                            femtobase::Px<femtobase::Pt, femtobase::Eta>,
                            femtobase::Py<femtobase::Pt, femtobase::Eta>,
                            femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                            femtobase::P<femtobase::Pt, femtobase::Eta>);
using FVzeros = FVzeros_001;

DECLARE_SOA_TABLE_VERSIONED(FVzeroMasks_001, "AOD", "FVZEROMASKS", 1,
                            femtovzeros::VzeroMask);
using FVzeroMasks = FVzeroMasks_001;

DECLARE_SOA_TABLE_VERSIONED(FVzeroExtras_001, "AOD", "FVZEROEXTRAS", 1,
                            femtovzeros::DauDCA,
                            femtovzeros::DecayVtxX,
                            femtovzeros::DecayVtxY,
                            femtovzeros::DecayVtxZ,
                            femtovzeros::TransRadius);
using FVzeroExtras = FVzeroExtras_001;

DECLARE_SOA_TABLE_VERSIONED(FVzeroPosDaus_001, "AOD", "FVZEROPOSDAUS", 1,
                            femtovzeros::TrackId,
                            femtobase::Pt,
                            femtobase::Eta,
                            femtobase::Phi,
                            femtovzeros::DauTrackMask,
                            femtovzeros::DauTPCMask);
using FVzeroPosDaus = FVzeroPosDaus_001;

DECLARE_SOA_TABLE_VERSIONED(FVzeroNegDaus_001, "AOD", "FVZERONEGDAUS", 1,
                            femtovzeros::TrackId,
                            femtobase::Pt,
                            femtobase::Eta,
                            femtobase::Phi,
                            femtovzeros::DauTrackMask,
                            femtovzeros::DauTPCMask);
using FVzeroNegDaus = FVzeroNegDaus_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOVZEROSDERIVED_H_
