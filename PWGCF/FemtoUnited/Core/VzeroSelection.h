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

/// \file FemtoDreamTrackCuts.h
/// \brief Definition of the FemtoDreamTrackCuts
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_

#include <cmath>
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited//Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace TrackSel
{
/// The different selections this task is capable of doing
enum TrackSel { kSign,         ///< Sign of the track
                kTPCnClsMin,   ///< Min. number of TPC clusters
                kTPCcRowsMin,  ///< Min. number of crossed TPC rows
                kTPCsClsMax,   ///< Max. number of shared TPC clusters
                kITSnClsMin,   ///< Min. number of ITS clusters
                kITSnClsIbMin, ///< Min. number of ITS clusters in the inner barrel
                kDCAxyMax,     ///< Max. DCA_xy (cm)
                kDCAzMax,      ///< Max. DCA_z (cm)
                kTrackselMax
};
} // namespace TrackSel

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TrackSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackMaskType, TrackSel::kTrackselMax>
{
 public:
  TrackSelection() {};
  virtual ~TrackSelection() = default;
  template <class track>
  void ApplySelections(track Track)
  {
    this->resetMinimalSelection();
    this->setBitmaskForObservable(TrackSel::kSign, Track.sign());
    this->setBitmaskForObservable(TrackSel::kTPCnClsMin, Track.tpcNClsFound());
    this->setBitmaskForObservable(TrackSel::kTPCcRowsMin, Track.tpcNClsCrossedRows());
    this->setBitmaskForObservable(TrackSel::kTPCsClsMax, Track.tpcNClsShared());
    this->setBitmaskForObservable(TrackSel::kITSnClsMin, Track.itsNCls());
    this->setBitmaskForObservable(TrackSel::kITSnClsIbMin, Track.itsNClsInnerBarrel());
    this->assembleBismask();
  };
}; // namespace femtoDream
} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_
