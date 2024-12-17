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

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKTOFSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKTOFSELECTION_H_

#include <cmath>
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited//Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace TrackTofSel
{
/// The different selections this task is capable of doing
enum TrackTpcSel { kElectron, ///< Electon PID
                   kPion,     ///< Pion PID
                   kKaon,     ///< Kaon PID
                   kProton,   ///< Proton PID
                   kDeuteron, ///< Deuteron PID
                   kTriton,   ///< Triton PID
                   kHelium3,  ///< He3 PID
                   kTrackTofSelMax
};
} // namespace TrackTofSel

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TrackTofSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackTOFMaskType, TrackTofSel::kTrackTofSelMax>
{
 public:
  TrackTofSelection() {};
  virtual ~TrackTofSelection() = default;
  template <class track>
  void ApplySelections(track Track)
  {
    this->setBitmaskForObservable(TrackTofSel::kElectron, Track.tofNSigmaEl());
    this->setBitmaskForObservable(TrackTofSel::kPion, Track.tofNSigmaPi());
    this->setBitmaskForObservable(TrackTofSel::kKaon, Track.tofNSigmaKa());
    this->setBitmaskForObservable(TrackTofSel::kProton, Track.tofNSigmaPr());
    this->setBitmaskForObservable(TrackTofSel::kDeuteron, Track.tofNSigmaDe());
    this->setBitmaskForObservable(TrackTofSel::kTriton, Track.tofNSigmaTr());
    this->setBitmaskForObservable(TrackTofSel::kHelium3, Track.tofNSigmaHe());
    this->assembleBismask();
  };
}; // namespace femtoDream
} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_
