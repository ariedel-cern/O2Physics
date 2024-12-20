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

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKTPCTOFSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKTPCTOFSELECTION_H_

#include <cmath>
#include "PWGCF/FemtoUnited/Utils/FemtoUtils.h"
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace tracktpctofselection
{
/// The different selections this task is capable of doing
enum TrackTpcTofSels { kElectron, ///< Electon PID
                       kPion,     ///< Pion PID
                       kKaon,     ///< Kaon PID
                       kProton,   ///< Proton PID
                       kDeuteron, ///< Deuteron PID
                       kTriton,   ///< Triton PID
                       kHelium,   ///< He PID
                       kTrackTpcTofSelMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TrackTpcTofSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackTPCTOFMaskType, TrackTpcTofSels::kTrackTpcTofSelMax>
{
 public:
  TrackTpcTofSelection() {}
  virtual ~TrackTpcTofSelection() = default;
  template <class track>
  void ApplySelections(track Track)
  {
    this->setBitmaskForObservable(TrackTpcTofSels::kElectron, utils::geometricMean(Track.tpcNSigmaEl(), Track.tofNSigmaEl()));
    this->setBitmaskForObservable(TrackTpcTofSels::kPion, utils::geometricMean(Track.tpcNSigmaPi(), Track.tofNSigmaPi()));
    this->setBitmaskForObservable(TrackTpcTofSels::kKaon, utils::geometricMean(Track.tpcNSigmaKa(), Track.tofNSigmaKa()));
    this->setBitmaskForObservable(TrackTpcTofSels::kProton, utils::geometricMean(Track.tpcNSigmaPr(), Track.tofNSigmaPr()));
    this->setBitmaskForObservable(TrackTpcTofSels::kDeuteron, utils::geometricMean(Track.tpcNSigmaDe(), Track.tofNSigmaDe()));
    this->setBitmaskForObservable(TrackTpcTofSels::kTriton, utils::geometricMean(Track.tpcNSigmaTr(), Track.tofNSigmaTr()));
    this->setBitmaskForObservable(TrackTpcTofSels::kHelium, utils::geometricMean(Track.tpcNSigmaHe(), Track.tofNSigmaHe()));
    this->assembleBismask();
  };
}; // namespace femtoDream
} // namespace tracktpctofselection
} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_TRACKTPCTOFSELECTION_H_
