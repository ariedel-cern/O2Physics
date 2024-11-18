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

/// \file CollisionHisto.h
/// \brief Femtounited collision histograms

#ifndef PWGCF_FEMTOUNITED_CORE_COLLISIONHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_COLLISIONHISTMANAGER_H_

#include <string_view>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <map>

#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace femtounitedcolhistmanager
{

enum EventVariable {
  kPosZ,
  kMult,
  kCent,
  kMagField,
  kPosX,
  kPosY,
  kEventVariableLast
};

constexpr std::string_view OutputDir = "CollisionHistograms/";

constexpr std::array<std::pair<std::string_view, std::string_view>, kEventVariableLast> HistogramNames = {
  {{"VertexZ", "Vertex Z; V_{Z} (cm); Entries"},
   {"Multiplicity", "Multiplicity; mult; Entries"},
   {"Centrality", "Centrality; cent (%); Entries"},
   {"MagneticField", "Magnetic Field; B (T); Entries"},
   {"VertexX", ""}}};

enum class Mode {
  kANALYSIS,
  kQA,
  kMC,
  kLAST_MODE
};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
class CollisionHistManager
{
 public:
  /// Destructor
  virtual ~CollisionHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  template <Mode mode>
  void init(HistogramRegistry* registry, std::map<EventVariable, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;
    std::string Dir = std::string(OutputDir);

    if constexpr (mode == Mode::kANALYSIS || mode == Mode::kQA) {
      mHistogramRegistry->add((Dir + std::string(HistogramNames[kPosZ].first)).c_str(), HistogramNames[kPosZ].second.data(), HistType::kTH1F, {Specs[kPosZ]});
      mHistogramRegistry->add((Dir + std::string(HistogramNames[kMult].first)).c_str(), HistogramNames[kMult].second.data(), HistType::kTH1F, {Specs[kMult]});
      mHistogramRegistry->add((Dir + std::string(HistogramNames[kCent].first)).c_str(), HistogramNames[kCent].second.data(), HistType::kTH1F, {Specs[kCent]});
      mHistogramRegistry->add((Dir + std::string(HistogramNames[kMagField].first)).c_str(), HistogramNames[kMagField].second.data(), HistType::kTH1F, {Specs[kMagField]});
    }
  }

  template <Mode mode, typename T>
  void fill(T const& col)
  {
    if constexpr (mode == Mode::kANALYSIS) {
      mHistogramRegistry->fill(HIST(OutputDir) + HIST(HistogramNames[kPosZ].first), col.posZ());
      mHistogramRegistry->fill(HIST(OutputDir) + HIST(HistogramNames[kMult].first), col.mult());
      mHistogramRegistry->fill(HIST(OutputDir) + HIST(HistogramNames[kCent].first), col.cent());
      mHistogramRegistry->fill(HIST(OutputDir) + HIST(HistogramNames[kMagField].first), col.magField());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;
}; // namespace femtounitedcolhistmanager
} // namespace femtounitedcolhistmanager

#endif // PWGCF_FEMTOUNITED_CORE_COLLISIONHISTMANAGER_H_
