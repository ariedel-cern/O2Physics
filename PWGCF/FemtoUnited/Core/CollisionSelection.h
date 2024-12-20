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

/// \file ColSelection.h
/// \brief Collision selection
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_COLSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_COLSELECTION_H_

#include <cmath>
#include "PWGCF/FemtoUnited/Core/Modes.h"

namespace o2::analysis::femtounited
{
namespace collisionselection
{
/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class CollisionSelection
{
 public:
  CollisionSelection() {}
  virtual ~CollisionSelection() = default;

  void init(bool OfflineSelection)
  {
    mOfflineSelection = OfflineSelection;
  };

  template <modes::System system, typename T>
  bool isSelected(T col)
  {
    bool flag = false;
    if constexpr (modes::isSystemSet(system, modes::System::kRun3)) {
      if (mOfflineSelection && col.sel8()) {
        flag = true;
      }
    }
    return flag;
  }

 private:
  bool mOfflineSelection = false;
};
}; // namespace collisionselection
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_COLSELECTION_H_
