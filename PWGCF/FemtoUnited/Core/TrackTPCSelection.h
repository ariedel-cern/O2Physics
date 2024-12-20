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
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKTPCSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKTPCSELECTION_H_

#include <cmath>
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited//Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace tracktpcselection
{
/// The different selections this task is capable of doing
enum TrackTpcSels { kElectron, ///< Electon PID
                    kPion,     ///< Pion PID
                    kKaon,     ///< Kaon PID
                    kProton,   ///< Proton PID
                    kDeuteron, ///< Deuteron PID
                    kTriton,   ///< Triton PID
                    kHelium,   ///< He3 PID
                    kTrackTpcSelMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TrackTpcSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackTPCMaskType, TrackTpcSels::kTrackTpcSelMax>
{
 public:
  TrackTpcSelection() {}
  virtual ~TrackTpcSelection() = default;
  template <class track>
  void ApplySelections(track Track)
  {
    this->setBitmaskForObservable(TrackTpcSels::kElectron, Track.tpcNSigmaEl());
    this->setBitmaskForObservable(TrackTpcSels::kPion, Track.tpcNSigmaPi());
    this->setBitmaskForObservable(TrackTpcSels::kKaon, Track.tpcNSigmaKa());
    this->setBitmaskForObservable(TrackTpcSels::kProton, Track.tpcNSigmaPr());
    this->setBitmaskForObservable(TrackTpcSels::kDeuteron, Track.tpcNSigmaDe());
    this->setBitmaskForObservable(TrackTpcSels::kTriton, Track.tpcNSigmaTr());
    this->setBitmaskForObservable(TrackTpcSels::kHelium, Track.tpcNSigmaHe());
    this->assembleBismask();
  };
}; // namespace femtoDream
}; // namespace tracktpcselection
}; // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_TRACKTPCSELECTION_H_
