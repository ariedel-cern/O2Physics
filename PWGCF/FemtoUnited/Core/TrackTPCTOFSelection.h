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
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited//Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace TrackTpcTofSel
{
/// The different selections this task is capable of doing
enum TrackTpcSel { kElectron, ///< Electon PID
                   kPion,     ///< Pion PID
                   kKaon,     ///< Kaon PID
                   kProton,   ///< Proton PID
                   kDeuteron, ///< Deuteron PID
                   kTriton,   ///< Triton PID
                   kHelium3,  ///< He3 PID
                   kTrackTpcTofSelMax
};
} // namespace TrackTpcTofSel

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TrackTpcTofSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackTPCTOFMaskType, TrackTpcTofSel::kTrackTpcTofSelMax>
{
 public:
  TrackTpcTofSelection() {};
  virtual ~TrackTpcTofSelection() = default;
  template <class track>
  void ApplySelections(track Track)
  {
    this->setBitmaskForObservable(TrackTpcTofSel::kElectron, std::sqrt(std::pow(Track.tpcNSigmaEl(), 2) + std::pow(Track.tofNSigmaEl(), 2)));
    this->setBitmaskForObservable(TrackTpcTofSel::kPion, std::sqrt(std::pow(Track.tpcNSigmaPi(), 2) + std::pow(Track.tofNSigmaPi(), 2)));
    this->setBitmaskForObservable(TrackTpcTofSel::kKaon, std::sqrt(std::pow(Track.tpcNSigmaKa(), 2) + std::pow(Track.tofNSigmaKa(), 2)));
    this->setBitmaskForObservable(TrackTpcTofSel::kProton, std::sqrt(std::pow(Track.tpcNSigmaPr(), 2) + std::pow(Track.tofNSigmaPr(), 2)));
    this->setBitmaskForObservable(TrackTpcTofSel::kDeuteron, std::sqrt(std::pow(Track.tpcNSigmaDe(), 2) + std::pow(Track.tofNSigmaDe(), 2)));
    this->setBitmaskForObservable(TrackTpcTofSel::kTriton, std::sqrt(std::pow(Track.tpcNSigmaTr(), 2) + std::pow(Track.tofNSigmaTr(), 2)));
    this->setBitmaskForObservable(TrackTpcTofSel::kHelium3, std::sqrt(std::pow(Track.tpcNSigmaHe(), 2) + std::pow(Track.tofNSigmaHe(), 2)));
    this->assembleBismask();
  };
}; // namespace femtoDream
} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_
