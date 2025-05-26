// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.femtobaseder
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoVzerosDerived.h
/// \brief v0 tables tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

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
DECLARE_SOA_COLUMN(VzeroMass, vzeroMass, float);         //! Mass of Vzero
DECLARE_SOA_COLUMN(AntiVzeroMass, antiVzeroMass, float); //! Mass of Vzero

// columns for Vzero bit masks
DECLARE_SOA_COLUMN(VzeroMask, vzeroMask, femtodatatypes::VzeroMaskType);                   //! Bitmask for Vzero selections
DECLARE_SOA_COLUMN(VzeroDauPidMask, vzeroDauPidMask, femtodatatypes::VzeroDauPidMaskType); //! Bitmask for Vzero selections

// columns for Vzero debug information
DECLARE_SOA_COLUMN(DauDCA, dauDCA, float);           //! Vzero daughter DCA at decay vertex
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Vzero transvers radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! x coordinate of Vzero decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! y coordinate of Vzero decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! z coordinate of Vzero decay vertex
DECLARE_SOA_COLUMN(KaonMass, kaonMass, float);       //! Vzero mass using Kaon hypothesis

// columns for Vzero daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int, FUTracks, "_Pos"); //!
DECLARE_SOA_COLUMN(PosDauPt, posDauPt, float);                        //! Positive daughter pt
DECLARE_SOA_COLUMN(PosDauEta, posDauEta, float);                      //! Positive daughter eta
DECLARE_SOA_COLUMN(PosDauPhi, posDauPhi, float);                      //! Positive daughter phi
DECLARE_SOA_DYNAMIC_COLUMN(PosDauP, posDauP,                          //! Momentum
                           [](float pt, float eta) -> float {
                             return pt * std::cosh(eta);
                           });

DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int, FUTracks, "_Neg"); //!
DECLARE_SOA_COLUMN(NegDauPt, negDauPt, float);                        //! Negative daughter pt
DECLARE_SOA_COLUMN(NegDauEta, negDauEta, float);                      //! Negative daughter eta
DECLARE_SOA_COLUMN(NegDauPhi, negDauPhi, float);                      //! Negative daughter phi
DECLARE_SOA_DYNAMIC_COLUMN(NegDauP, negDauP,                          //! Momentum
                           [](float pt, float eta) -> float {
                             return pt * std::cosh(eta);
                           });

// columns for Vzero daughter debug information
DECLARE_SOA_COLUMN(PosDauTpcNclsFound, posDauTpcNclsFound, uint8_t);     //!
DECLARE_SOA_COLUMN(PosDauDcaxy, posDauDcaxy, float);                     //!
DECLARE_SOA_COLUMN(PosDauDcaz, posDauDcaz, float);                       //!
DECLARE_SOA_COLUMN(PosDauProtonTpcNsigma, posDauProtonTpcNsigma, float); //!
DECLARE_SOA_COLUMN(PosDauPionTpcNsigma, posDauPionTpcNsigma, float);     //!
DECLARE_SOA_COLUMN(NegDauTpcNclsFound, negDauTpcNclsFound, uint8_t);     //!
DECLARE_SOA_COLUMN(NegDauDcaxy, negDauDcaxy, float);                     //!
DECLARE_SOA_COLUMN(NegDauDcaz, negDauDcaz, float);                       //!
DECLARE_SOA_COLUMN(NegDauProtonTpcNsigma, negDauProtonTpcNsigma, float); //!
DECLARE_SOA_COLUMN(NegDauPionTpcNsigma, negDauPionTpcNsigma, float);     //!

} // namespace femtovzeros

// table for basic vzero information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUVzeros_001, "FUVZEROS", 1,
                                   o2::soa::Index<>,
                                   femtobase::CollisionId,
                                   femtobase::Pt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtovzeros::VzeroMass,
                                   femtovzeros::AntiVzeroMass,
                                   femtobase::Theta<femtobase::Eta>,
                                   femtobase::Px<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Py<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                                   femtobase::P<femtobase::Pt, femtobase::Eta>);
using FUVzeros = FUVzeros_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUVzeroMasks_001, "FUVZEROMASKS", 1,
                                   femtovzeros::VzeroMask,
                                   femtovzeros::VzeroDauPidMask);
using FUVzeroMasks = FUVzeroMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUVzeroExtras_001, "FUVZEROEXTRAS", 1,
                                   femtovzeros::DauDCA,
                                   femtovzeros::DecayVtxX,
                                   femtovzeros::DecayVtxY,
                                   femtovzeros::DecayVtxZ,
                                   femtovzeros::TransRadius,
                                   femtovzeros::KaonMass);
using FUVzeroExtras = FUVzeroExtras_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUVzeroDaus_001, "FUVZERODAUS", 1,
                                   femtovzeros::PosDauId,
                                   femtovzeros::PosDauPt,
                                   femtovzeros::PosDauEta,
                                   femtovzeros::PosDauPhi,
                                   femtovzeros::NegDauId,
                                   femtovzeros::NegDauPt,
                                   femtovzeros::NegDauEta,
                                   femtovzeros::NegDauPhi,
                                   femtovzeros::PosDauP<femtovzeros::PosDauPt, femtovzeros::PosDauEta>,
                                   femtovzeros::NegDauP<femtovzeros::NegDauPt, femtovzeros::NegDauEta>);
using FUVzeroDaus = FUVzeroDaus_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUVzeroDauExts_001, "FUVZERODAUEXTS", 1,
                                   femtovzeros::PosDauTpcNclsFound,
                                   femtovzeros::PosDauDcaxy,
                                   femtovzeros::PosDauDcaz,
                                   femtovzeros::PosDauProtonTpcNsigma,
                                   femtovzeros::PosDauPionTpcNsigma,
                                   femtovzeros::NegDauTpcNclsFound,
                                   femtovzeros::NegDauDcaxy,
                                   femtovzeros::NegDauDcaz,
                                   femtovzeros::NegDauProtonTpcNsigma,
                                   femtovzeros::NegDauPionTpcNsigma);
using FUVzeroDauExts = FUVzeroDauExts_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOVZEROSDERIVED_H_
