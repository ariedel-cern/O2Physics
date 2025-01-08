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

/// \file CollisionSelection.h
/// \brief collision selection
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_COLLISIONSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_COLLISIONSELECTION_H_

#include <cmath>
#include <string>
#include "Framework/Configurable.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

namespace o2::analysis::femtounited
{
namespace collisionselection
{

struct ConfCollisionSelection : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("CollisionSelection");
  o2::framework::Configurable<float> vtxZMin{"vtxZMin", -10., "Minimum vertex Z position (cm)"};
  o2::framework::Configurable<float> vtxZMax{"vtxZMax", 10., "Maximum vertex Z position (cm)"};
  o2::framework::Configurable<float> multMin{"multMin", 0, "Minimum multiplicity"};
  o2::framework::Configurable<float> multMax{"multMax", 200, "Maximum multiplicity"};
  o2::framework::Configurable<float> centMin{"centMin", 0.0f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> centMax{"centMax", 100.0f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> spherMin{"spherMin", 0.0f, "Minimum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> spherMax{"spherMax", 2.0f, "Maximum centrality (multiplicity percentile)"};
  o2::framework::Configurable<float> magFieldMin{"magFieldMin", -1.0f, "Minimum magnetic field strength (T)"};
  o2::framework::Configurable<float> magFieldMax{"magFieldMax", 1.0f, "Maximum magnetic field strength (T)"};
};

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
  bool isSelected(T const& col)
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
#endif // PWGCF_FEMTOUNITED_CORE_COLLISIONSELECTION_H_
