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

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTRACKSDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTRACKSDERIVED_H_

#include <cmath>
#include "Framework/ASoA.h"
#include "Framework/Expressions.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"

namespace o2::aod
{
namespace femtotracks
{

// columns for track selections
DECLARE_SOA_COLUMN(TrackMask, trackMask, femtodatatypes::TrackMaskType);         //! Bitmask for track selections
DECLARE_SOA_COLUMN(TpcMask, tpcMask, femtodatatypes::TrackTPCMaskType);          //! Bitmask for TPC PID selections
DECLARE_SOA_COLUMN(TofMask, tofMask, femtodatatypes::TrackTOFMaskType);          //! Bitmask for TOF PID selections
DECLARE_SOA_COLUMN(TpcTofMask, tpctofMask, femtodatatypes::TrackTPCTOFMaskType); //! Bitmask for combined TPC+TOF PID selections

// columns for DCA
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);                                                                                       //! Dca in XY plane
DECLARE_SOA_COLUMN(DcaZ, dcaZ, float);                                                                                         //! Dca in Z direction
DECLARE_SOA_DYNAMIC_COLUMN(Dca, dca, [](float dcaXY, float dcaZ) -> float { return std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ); }); //! Dca

// columns for track debug information
DECLARE_SOA_COLUMN(Sign, sign, int8_t); //! Sign (charge)

// its related information
DECLARE_SOA_COLUMN(IsPVContributor, isPVContributor, bool);          //! True if track is PV contributer
DECLARE_SOA_COLUMN(ITSNCls, itsNCls, uint8_t);                       //! Number of ITS clusters (max 7)
DECLARE_SOA_COLUMN(ITSNClsInnerBarrel, itsNClsInnerBarrel, uint8_t); //! Number of ITS clusters in the inner barrel (max 3)
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, float);                   //! ITS chi2 / cluster
DECLARE_SOA_COLUMN(ITSClusterSizes, itsclusterSizes, uint32_t);      //! ITS cluster sizes (4 bits per layer)

// tpc related information
DECLARE_SOA_COLUMN(TPCSignal, tpcSignal, float);                                   //! TPC signal
DECLARE_SOA_COLUMN(TPCInnerParam, tpcInnerParam, bool);                            //! Momentum at inner wall of TPC
DECLARE_SOA_COLUMN(TPCNClsFindable, tpcNClsFindable, uint8_t);                     //! Number of TPC clusters
DECLARE_SOA_COLUMN(TPCNClsCrossedRows, tpcNClsCrossedRows, uint8_t);               //! Number of TPC crossed rows
DECLARE_SOA_DYNAMIC_COLUMN(TPCCrossedRowsOverFindable, tpcCrossedRowsOverFindable, //! Number of crossed rows over found TPC clusters
                           [](uint8_t tpcNClsFindable, uint8_t tpcNClsCrossedRows) -> float { return (float)tpcNClsCrossedRows / (float)tpcNClsFindable; });
DECLARE_SOA_COLUMN(TPCNClsShared, tpcNClsShared, uint8_t); //! Number of shared TPC clusters
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, float);         //! TPC chi2 / findable clusters

// tof related information
DECLARE_SOA_COLUMN(TOFBeta, tofbeta, float); //! TOF beta

// TPC PID information
DECLARE_SOA_COLUMN(TPCNSigmaEl, tpcNSigmaEl, float);   //! Nsigma separation with the TPC detector for electron
DECLARE_SOA_COLUMN(TPCNSigmaPi, tpcNSigmaPi, float);   //! Nsigma separation with the TPC detector for pion
DECLARE_SOA_COLUMN(TPCNSigmaKa, tpcNSigmaKa, float);   //! Nsigma separation with the TPC detector for kaon
DECLARE_SOA_COLUMN(TPCNSigmaPr, tpcNSigmaPr, float);   //! Nsigma separation with the TPC detector for proton
DECLARE_SOA_COLUMN(TPCNSigmaDe, tpcNSigmaDe, float);   //! Nsigma separation with the TPC detector for deuteron
DECLARE_SOA_COLUMN(TPCNSigmaTr, tpcNSigmaTr, float);   //! Nsigma separation with the TPC detector for triton
DECLARE_SOA_COLUMN(TPCNSigmaHe3, tpcNSigmaHe3, float); //! Nsigma separation with the TPC detector for helium3

// TOF PID information
DECLARE_SOA_COLUMN(TOFNSigmaEl, tofNSigmaEl, float);   //! Nsigma separation with the TOF detector for electron
DECLARE_SOA_COLUMN(TOFNSigmaPi, tofNSigmaPi, float);   //! Nsigma separation with the TOF detector for pion
DECLARE_SOA_COLUMN(TOFNSigmaKa, tofNSigmaKa, float);   //! Nsigma separation with the TOF detector for kaon
DECLARE_SOA_COLUMN(TOFNSigmaPr, tofNSigmaPr, float);   //! Nsigma separation with the TOF detector for proton
DECLARE_SOA_COLUMN(TOFNSigmaDe, tofNSigmaDe, float);   //! Nsigma separation with the TOF detector for deuteron
DECLARE_SOA_COLUMN(TOFNSigmaTr, tofNSigmaTr, float);   //! Nsigma separation with the TOF detector for triton
DECLARE_SOA_COLUMN(TOFNSigmaHe3, tofNSigmaHe3, float); //! Nsigma separation with the TOF detector for helium3
} // namespace femtotracks

// table for basic track information
DECLARE_SOA_TABLE_VERSIONED(FUTracks_001, "AOD", "FUTRACKS", 1,
                            o2::soa::Index<>,
                            femtobase::CollisionId,
                            femtobase::Pt,
                            femtobase::Eta,
                            femtobase::Phi,
                            femtobase::Theta<femtobase::Eta>,
                            femtobase::Px<femtobase::Pt, femtobase::Eta>,
                            femtobase::Py<femtobase::Pt, femtobase::Eta>,
                            femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                            femtobase::P<femtobase::Pt, femtobase::Eta>);
using FUTracks = FUTracks_001;

// table for track selections and PID selections
DECLARE_SOA_TABLE_VERSIONED(FUTrackMasks_001, "AOD", "FUTRACKMASKS", 1,
                            femtotracks::TrackMask,
                            femtotracks::TpcMask,
                            femtotracks::TofMask,
                            femtotracks::TpcTofMask);
using FUTrackMasks = FUTrackMasks_001;

// table for track DCA
DECLARE_SOA_TABLE_VERSIONED(FUTrackDCAs_001, "AOD", "FUTRACKDCAS", 1,
                            femtotracks::DcaXY,
                            femtotracks::DcaZ,
                            femtotracks::Dca<femtotracks::DcaXY, femtotracks::DcaZ>);
using FUTrackDCAs = FUTrackDCAs_001;

// table for extra track information
DECLARE_SOA_TABLE_VERSIONED(FUTrackExtras_001, "AOD", "FUTRACKEXTRAS", 1,
                            femtotracks::Sign,
                            femtotracks::IsPVContributor,
                            femtotracks::ITSNCls,
                            femtotracks::ITSNClsInnerBarrel,
                            femtotracks::ITSChi2NCl,
                            femtotracks::ITSClusterSizes,
                            femtotracks::TPCSignal,
                            femtotracks::TPCInnerParam,
                            femtotracks::TPCNClsFindable,
                            femtotracks::TPCNClsCrossedRows,
                            femtotracks::TPCNClsShared,
                            femtotracks::TOFBeta,
                            femtotracks::TPCCrossedRowsOverFindable<femtotracks::TPCNClsFindable, femtotracks::TPCNClsCrossedRows>);
using FUTrackExtras = FUTrackExtras_001;

// table for extra PID information
DECLARE_SOA_TABLE_VERSIONED(FUTrackPids_001, "AOD", "FUTRACKPIDS", 1,
                            femtotracks::TPCNSigmaEl,
                            femtotracks::TPCNSigmaPi,
                            femtotracks::TPCNSigmaKa,
                            femtotracks::TPCNSigmaPr,
                            femtotracks::TPCNSigmaDe,
                            femtotracks::TPCNSigmaTr,
                            femtotracks::TPCNSigmaHe3,
                            femtotracks::TOFNSigmaEl,
                            femtotracks::TOFNSigmaPi,
                            femtotracks::TOFNSigmaKa,
                            femtotracks::TOFNSigmaPr,
                            femtotracks::TOFNSigmaDe,
                            femtotracks::TOFNSigmaTr,
                            femtotracks::TOFNSigmaHe3);
using FUTrackPids = FUTrackPids_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTRACKSDERIVED_H_
