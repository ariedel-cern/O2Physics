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

/// \file VzeroSelection.h
/// \brief Vzero selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_

#include <algorithm>
#include <cmath>

#include "Framework/Configurable.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/BaseSelection.h"

namespace o2::analysis::femtounited
{
namespace vzeroselection
{
/// The different selections this task is capable of doing
enum VzeroSels {
  kDcaDaughMax,     ///< Max. DCA of the daughers at decay vertex
  kCpaMin,          ///< Min. CPA (cosine pointing angle)
  kTransRadMin,     ///< Min. transverse radius
  kTransRadMax,     ///< max. transverse radius
  kDecayVtxMax,     ///< Max. distance of decay vertex in x,y,z
  kPosDauDcaMin,    ///< Min. DCA of the positive daughers at primary vertex
  kPosDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter
  kNegDauDcaMin,    ///< Min. DCA of the positive daughers at primary vertex
  kNegDauTpcClsMin, ///< Min. number of TPC clusters of positive daughter
  kVzeroSelsMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class VzeroSelection : public BaseSelection<float, o2::aod::femtodatatypes::VzeroMaskType, VzeroSels::kVzeroSelsMax>
{
 public:
  VzeroSelection() {}
  virtual ~VzeroSelection() = default;

  void setKaonMassLimits(float lower, float upper)
  {
    if (lower < 0 || upper < 0) {
      mMassKaonLowerLimit = -1.f;
      mMassKaonUpperLimit = -1.f;
      return;
    }
    mMassKaonLowerLimit = lower;
    mMassKaonUpperLimit = upper;
  }

  template <typename T>
  bool checkKaonMassLimit(T const& v0)
  {
    if (mMassKaonLowerLimit < 0 || mMassKaonUpperLimit < 0) {
      return true;
    }
    if (v0.mK0Short() > mMassKaonLowerLimit && v0.mK0Short() < mMassKaonUpperLimit) {
      return false;
    } else {
      return true;
    }
  }

  template <class V0, class Tracks>
  void applySelections(V0 const& v0, Tracks const& /*tracks*/)
  {
    this->reset();
    // vzero selections
    this->setBitmaskForObservable(VzeroSels::kDcaDaughMax, v0.dcaV0daughters());
    this->setBitmaskForObservable(VzeroSels::kCpaMin, v0.v0cosPA());
    this->setBitmaskForObservable(VzeroSels::kTransRadMin, v0.v0radius());
    this->setBitmaskForObservable(VzeroSels::kTransRadMax, v0.v0radius());
    // for decay vertex, the x,y and z coordinate have to be below a certain threshold
    // compare the largest of the 3 to the limit to set the bit
    std::array<float, 3> decayCoordinates = {v0.x(), v0.y(), v0.z()};
    this->setBitmaskForObservable(VzeroSels::kDecayVtxMax, *std::max_element(decayCoordinates.begin(), decayCoordinates.end()));

    // positive daughter selections
    auto posDaughter = v0.template posTrack_as<Tracks>();
    this->setBitmaskForObservable(VzeroSels::kPosDauDcaMin, std::hypot(posDaughter.dcaXY(), posDaughter.dcaZ()));
    this->setBitmaskForObservable(VzeroSels::kPosDauTpcClsMin, posDaughter.tpcNClsFound());

    // negative daughter selections
    auto negDaughter = v0.template negTrack_as<Tracks>();
    this->setBitmaskForObservable(VzeroSels::kNegDauDcaMin, std::hypot(negDaughter.dcaXY(), negDaughter.dcaZ()));
    this->setBitmaskForObservable(VzeroSels::kNegDauTpcClsMin, negDaughter.tpcNClsFound());

    this->assembleBismask();
  };

 protected:
  float mMassKaonLowerLimit = -1.f;
  float mMassKaonUpperLimit = -1.f;
};
} // namespace vzeroselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_
