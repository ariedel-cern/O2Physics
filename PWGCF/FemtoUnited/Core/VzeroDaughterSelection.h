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
  kDcaMin,      ///< Min. DCA of the daughers at primary vertex
  kTpcClsMin,   ///< Min. transverse radius
  kDecayVtxMax, ///< Max. distance of decay vertex
  kVzeroSelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class VzeroDaughterSelection : public BaseSelection<float, o2::aod::femtodatatypes::VzeroDauTrackMaskType, VzeroDaughterSels::kVzeroSelsMax>
{
 public:
  VzeroDaughterSelection() {}
  virtual ~VzeroDaughterSelection() = default;
  // template <class V0>
  // void ApplySelections(V0 v0) {
  // this->resetMinimalSelection();
  // this->setBitmaskForObservable(VzeroSels::kDcaDaughMax, v0.dcaV0daughters());
  // this->setBitmaskForObservable(VzeroSels::kCpaMin, v0.v0cosPA());
  // this->setBitmaskForObservable(VzeroSels::kTransRadMin, v0.v0radius());
  // this->setBitmaskForObservable(VzeroSels::kTransRadMax, v0.v0radius());
  // this->assembleBismask();
  // };
};
} // namespace vzerodaughterselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_
