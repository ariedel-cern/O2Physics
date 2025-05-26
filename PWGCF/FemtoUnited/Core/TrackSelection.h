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

/// \file TrackSelection.h
/// \brief Definition of track selections
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_

#include <cmath>
#include "Framework/Configurable.h"
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace trackselection
{
/// The different selections this task is capable of doing
enum TrackSels {
  kTPCnClsMin,   ///< Min. number of TPC clusters
  kTPCcRowsMin,  ///< Min. number of crossed TPC rows
  kTPCsClsMax,   ///< Max. number of shared TPC clusters
  kITSnClsMin,   ///< Min. number of ITS clusters
  kITSnClsIbMin, ///< Min. number of ITS clusters in the inner barrel
  kDCAxyMax,     ///< Max. DCA_xy (cm) as a function of pT
  kDCAzMax,      ///< Max. DCA_z (cm) as a function of pT
  kTrackselMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TrackSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackMaskType, TrackSels::kTrackselMax>
{
 public:
  TrackSelection() {}
  virtual ~TrackSelection() = default;
  template <class track>
  void applySelections(track const& Track)
  {
    this->reset();
    this->setBitmaskForObservable(TrackSels::kTPCnClsMin, Track.tpcNClsFound());
    this->setBitmaskForObservable(TrackSels::kTPCcRowsMin, Track.tpcNClsCrossedRows());
    this->setBitmaskForObservable(TrackSels::kTPCsClsMax, Track.tpcNClsShared());
    this->setBitmaskForObservable(TrackSels::kITSnClsMin, Track.itsNCls());
    this->setBitmaskForObservable(TrackSels::kITSnClsIbMin, Track.itsNClsInnerBarrel());

    // evalue bitmask for pt dependent dca cuts
    this->updateLimits(TrackSels::kDCAxyMax, Track.pt());
    this->setBitmaskForObservable(TrackSels::kDCAxyMax, Track.dcaXY());

    this->updateLimits(TrackSels::kDCAzMax, Track.pt());
    this->setBitmaskForObservable(TrackSels::kDCAzMax, Track.dcaZ());

    this->assembleBismask();
  };
};
}; // namespace trackselection
}; // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_
