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

/// \file VzeroDaughterPidSelection.h
/// \brief track pid selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_VZERODAUGHTERPIDSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_VZERODAUGHTERPIDSELECTION_H_

#include <cmath>

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

namespace o2::analysis::femtounited
{
namespace vzerodaughterpidselection
{
/// The different selections this task is capable of doing
enum VzeroDaughterPidSels {

  // charge combination for lamdab
  kNegDaughTpcPion,   ///< TPC Pion PID for negative daughter
  kPosDaughTpcProton, ///< TPC Proton PID for positive daughter

  // charge combination for antilambda
  kPosDaughTpcPion,   ///< TPC Pion PID for positive daughter
  kNegDaughTpcProton, ///< TPC Proton PID for negative daughter

  VzeroDaughterPidSelsMax
};

/// \class VzeroDaughterPidSelectionPos
/// \brief Cut class to contain and execute all cuts applied to pid of vzerodaughters
class VzeroDaughterPidSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackPidMaskType, VzeroDaughterPidSels::VzeroDaughterPidSelsMax>
{
 public:
  VzeroDaughterPidSelection() {}
  virtual ~VzeroDaughterPidSelection() = default;
  template <class track>
  void applySelections(track const& posDaughter, track const& negDaughter)
  {
    this->reset();

    // apply bits for lambda
    this->setBitmaskForObservable(vzerodaughterpidselection::kNegDaughTpcPion, negDaughter.tpcNSigmaPi());
    this->setBitmaskForObservable(vzerodaughterpidselection::kPosDaughTpcProton, posDaughter.tpcNSigmaPr());

    // apply bits for antilambda
    this->setBitmaskForObservable(vzerodaughterpidselection::kPosDaughTpcPion, posDaughter.tpcNSigmaPi());
    this->setBitmaskForObservable(vzerodaughterpidselection::kNegDaughTpcProton, negDaughter.tpcNSigmaPr());

    this->assembleBismask();
  };

  bool isLambda()
  {
    return this->getAnySelection(vzerodaughterpidselection::kNegDaughTpcPion) && this->getAnySelection(vzerodaughterpidselection::kPosDaughTpcProton);
  }

  bool isAntiLambda()
  {
    return this->getAnySelection(vzerodaughterpidselection::kPosDaughTpcPion) && this->getAnySelection(vzerodaughterpidselection::kNegDaughTpcPion);
  }

}; // namespace femtoDream
}; // namespace vzerodaughterpidselection
}; // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_VZERODAUGHTERPIDSELECTION_H_
