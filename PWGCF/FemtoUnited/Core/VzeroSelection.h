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
#include "PWGCF/FemtoUnited/Utils/FemtoUtils.h"

namespace o2::analysis::femtounited
{
namespace vzeroselection
{
/// The different selections this task is capable of doing
enum VzeroSels {
  kDcaDaughMax,        ///< Max. DCA of the daughers at decay vertex
  kCpaMin,             ///< Min. CPA (cosine pointing angle)
  kTransRadMin,        ///< Min. transverse radius
  kTransRadMax,        ///< max. transverse radius
  kDecayVtxMax,        ///< Max. distance of decay vertex in x,y,z
  kPosDauDcaMin,       ///< Min. DCA of the positive daughers at primary vertex
  kPosDauTpcClsMin,    ///< Min. number of TPC clusters of positive daughter
  kPosDauTpcNsigmaMax, ///< Max. |nsigma| TPC of the positive daughter
  kNegDauDcaMin,       ///< Min. DCA of the positive daughers at primary vertex
  kNegDauTpcClsMin,    ///< Min. number of TPC clusters of positive daughter
  kNegDauTpcNsigmaMax, ///< Max. |nsigma| TPC of the positive daughter
  kSign,               ///< +1 for particle, -1 for antiparticle
                       /// IMPORTANT: Always let the sign be the last bit
                       /// Like this the sign will be added last to the bit and ends up on the least significant bits this means we can always select v0/antiv0 very easily
                       ///  XXXXXX01 <- antivzero
                       ///  XXXXXX10 <- v0
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

  float getLoosestPosDauPidSelection()
  {
    return mSelections.at(kPosDauTpcNsigmaMax).getLoosestSelection();
  }

  float getLoosestNegDauPidSelection()
  {
    return mSelections.at(kNegDauTpcNsigmaMax).getLoosestSelection();
  }

  template <class V0, class Tracks>
  void setV0Type(V0 const& v0, Tracks const& /*track*/)
  {
    // get daughter tracks
    auto posDaughter = v0.template posTrack_as<Tracks>();
    auto negDaughter = v0.template negTrack_as<Tracks>();

    // get loosest pid selection
    float pidPosDaughter = getLoosestPosDauPidSelection();
    float pidNegDaughter = getLoosestNegDauPidSelection();

    // check with pid selection if v0 is lambda or antilambda
    bool isLambda = false;
    bool isAntiLambda = false;
    if (std::abs(posDaughter.tpcNSigmaPr()) <= pidPosDaughter && std::abs(negDaughter.tpcNSigmaPi()) < pidNegDaughter) {
      isLambda = true;
    }
    if (std::abs(posDaughter.tpcNSigmaPi()) <= pidPosDaughter && std::abs(negDaughter.tpcNSigmaPr()) < pidNegDaughter) {
      isAntiLambda = true;
    }

    // if both are false, return 0
    // this should cause the minimal selection to be false since the sign should be check to be either +1 or -1
    if (isLambda == false && isAntiLambda == false) {
      return;
    }

    // check that only one options was selected, if not we need to check further
    if (isLambda == true && isAntiLambda == false) {
      mSign = 1;
      mMass = v0.mLambda();
      return;
    }
    if (isLambda == false && isAntiLambda == true) {
      mSign = -1;
      mMass = v0.mAntiLambda();
      return;
    }

    // if PID is not enough, then we take the one closest to nominal mass
    float diffLambda = std::abs(o2::constants::physics::MassLambda - v0.mLambda());
    float diffAntiLambda = std::abs(o2::constants::physics::MassLambda - v0.mAntiLambda());

    if (diffLambda <= diffAntiLambda) {
      mMass = v0.mLambda();
      mSign = 1;
      return;
    } else {
      mMass = v0.mAntiLambda();
      mSign = -1;
      return;
    }
  };

  int getSign() { return mSign; }
  float getMass() { return mMass; }
  float getPosDaughterTpcNsigma() { return mPosDaughterTpcNsigma; }
  float getNegDaughterTpcNsigma() { return mNegDaughterTpcNsigma; }

  template <class V0, class Tracks>
  void applySelections(V0 const& v0, Tracks const& tracks)
  {
    // reset variables
    mMass = 0.f;
    mSign = 0;
    mPosDaughterTpcNsigma = -99.f;
    mNegDaughterTpcNsigma = -99.f;

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
    this->setBitmaskForObservable(VzeroSels::kPosDauDcaMin, utils::dca(posDaughter.dcaXY(), posDaughter.dcaZ()));
    this->setBitmaskForObservable(VzeroSels::kPosDauTpcClsMin, posDaughter.tpcNClsFound());

    // positive daughter selections
    auto negDaughter = v0.template negTrack_as<Tracks>();
    this->setBitmaskForObservable(VzeroSels::kNegDauDcaMin, utils::dca(negDaughter.dcaXY(), negDaughter.dcaZ()));
    this->setBitmaskForObservable(VzeroSels::kNegDauTpcClsMin, negDaughter.tpcNClsFound());

    // now we need to figure out if the v0 is a lambda or antilambda
    // this information is stored in the sign (similar to tracks)
    // +1 for particle and -1 for antiparticle
    this->setV0Type(v0, tracks);
    this->setBitmaskForObservable(VzeroSels::kSign, mSign);
    if (mSign == 1) {
      mPosDaughterTpcNsigma = posDaughter.tpcNSigmaPr();
      mNegDaughterTpcNsigma = negDaughter.tpcNSigmaPi();
    }
    if (mSign == -1) {
      mPosDaughterTpcNsigma = posDaughter.tpcNSigmaPi();
      mNegDaughterTpcNsigma = negDaughter.tpcNSigmaPr();
    }

    this->setBitmaskForObservable(VzeroSels::kPosDauTpcNsigmaMax, mPosDaughterTpcNsigma);
    this->setBitmaskForObservable(VzeroSels::kPosDauTpcNsigmaMax, mNegDaughterTpcNsigma);

    this->assembleBismask();
  };

 protected:
  int mSign = 0;
  float mMass = 0.f;
  float mMassKaonLowerLimit = -1.f;
  float mMassKaonUpperLimit = -1.f;
  float mPosDaughterTpcNsigma = -99.f;
  float mNegDaughterTpcNsigma = -99.f;
};
} // namespace vzeroselection
} // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_VZEROSELECTION_H_
