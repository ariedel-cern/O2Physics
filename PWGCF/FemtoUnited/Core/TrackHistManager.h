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

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKHISTMANAGER_H_

#include <string_view>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <map>

#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace femtounitedtrackhistmanager
{

enum TrackVariable {
  kPt,
  kEta,
  kPhi,
  kEventVariableLast
};

constexpr std::string_view OutputDir = "TrackHistograms/";

constexpr std::array<std::pair<std::string_view, std::string_view>, kEventVariableLast> HistogramNames = {
  {{"hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {"hEta", "Pseudorapdity; #eta; Entries"},
   {"hPhi", "Azimuthal angle; #varphi; Entries"}}};

enum class Mode {
  kANALYSIS,
  kQA,
  kMC,
  kLAST_MODE
};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
class TrackHistManager
{
 public:
  /// Destructor
  virtual ~TrackHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  template <Mode mode>
  void init(HistogramRegistry* registry, std::map<TrackVariable, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;
    std::string Dir = std::string(OutputDir);

    if constexpr (mode == Mode::kANALYSIS || mode == Mode::kQA) {
      mHistogramRegistry->add((Dir + std::string(HistogramNames[kPt].first)).c_str(), HistogramNames[kPt].second.data(), HistType::kTH1F, {Specs[kPt]});
      mHistogramRegistry->add((Dir + std::string(HistogramNames[kEta].first)).c_str(), HistogramNames[kEta].second.data(), HistType::kTH1F, {Specs[kEta]});
      mHistogramRegistry->add((Dir + std::string(HistogramNames[kPhi].first)).c_str(), HistogramNames[kPhi].second.data(), HistType::kTH1F, {Specs[kPhi]});
    }
  }

  template <Mode mode, typename T>
  void fill(T const& track)
  {
    if constexpr (mode == Mode::kANALYSIS) {
      mHistogramRegistry->fill(HIST(OutputDir) + HIST(HistogramNames[kPt].first), track.pt());
      mHistogramRegistry->fill(HIST(OutputDir) + HIST(HistogramNames[kEta].first), track.eta());
      mHistogramRegistry->fill(HIST(OutputDir) + HIST(HistogramNames[kPhi].first), track.phi());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;
}; // namespace femtounitedcolhistmanager
} // namespace femtounitedtrackhistmanager

#endif // PWGCF_FEMTOUNITED_CORE_COLLISIONHISTMANAGER_H_
