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

#ifndef PWGCF_FEMTOUNITED_CORE_VZERODAUGHTERSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_VZERODAUGHTERSELECTION_H_

#include <cmath>
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited//Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace vzerodaughterselection
{
/// The different selections this task is capable of doing
enum VzeroDaughterSels {
  kDcaMin,    ///< Min. DCA of the daughers at primary vertex
  kTpcClsMin, ///< Min. transverse radius
  kPid,       ///< Max. |nsigma tpc| for proton/pion pid for positive/negative daugher of V0/anti-V0
  kVzerosDaughterSelMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class VzeroDaughterSelection : public BaseSelection<float, o2::aod::femtodatatypes::VzeroDaughterMaskType, VzeroDaughterSels::kVzerosDaughterSelMax>
{
 public:
  VzeroDaughterSelection() {}
  virtual ~VzeroDaughterSelection() = default;
  template <class track>
  void ApplySelections(track const& Track, int signVzero, int signDaughter)
  {
    this->resetMinimalSelection();
    this->setBitmaskForObservable(VzeroDaughterSels::kDcaMin, Track.dcaXY());
    this->setBitmaskForObservable(VzeroDaughterSels::kTpcClsMin, Track.tpcNClsFound());

    // PID for daugher
    switch (signVzero) {
      case 1:
        if (signDaughter == 1) {
          this->setBitmaskForObservable(VzeroDaughterSels::kPid, Track.tpcNSigmaPr());
        }
        if (signDaughter == -1) {
          this->setBitmaskForObservable(VzeroDaughterSels::kPid, Track.tpcNSigmaPi());
        }
        break;
      case -1:
        // PID for daughter of anti-V0
        if (signDaughter == -1) {
          this->setBitmaskForObservable(VzeroDaughterSels::kPid, Track.tpcNSigmaPr());
        }
        if (signDaughter == 1) {
          this->setBitmaskForObservable(VzeroDaughterSels::kPid, Track.tpcNSigmaPi());
        }
        break;
      default:
        this->setBitmaskForObservable(VzeroDaughterSels::kPid, -999);
    }

    this->assembleBismask();
  };
};
} // namespace vzerodaughterselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_
