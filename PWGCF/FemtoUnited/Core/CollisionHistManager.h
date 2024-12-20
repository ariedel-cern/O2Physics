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

/// \file ColHistManager.h
/// \brief Femtounited collision histograms

#ifndef PWGCF_FEMTOUNITED_CORE_COLLISIONHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_COLLISIONHISTMANAGER_H_

#include <string_view>
#include <string>
#include <vector>
#include <array>
#include <map>

#include "Framework/HistogramRegistry.h"

#include "PWGCF/FemtoUnited/Core/HistManager.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

using namespace o2::framework;

namespace o2::analysis::femtounited
{
namespace colhistmanager
{

enum ColHist {
  kPosz,
  kMult,
  kCent,
  kMagField,
  kSphericity,
  // 2d qa
  kPoszVsMult,
  kPoszVsCent,
  kCentVsMult,
  kCentVsSphericity,
  kMultVsSphericity,
  kColHistLast
};

constexpr std::string_view ColAnalysisDir = "CollisionHistograms/Analysis/";
constexpr std::string_view ColQaDir = "CollisionHistograms/QA/";

constexpr std::array<Histmanager::HistInfo<ColHist>, kColHistLast> HistTable = {
  {
    {kPosz, kTH1F, "hPosz", "Vertex Z; V_{Z} (cm); Entries"},
    {kMult, kTH1F, "hMult", "Multiplicity; Multiplicity; Entries"},
    {kCent, kTH1F, "hCent", "Centrality; Centrality (%); Entries"},
    {kMagField, kTH1F, "hMagField", "Magnetic Field; B (T); Entries"},
    {kSphericity, kTH1F, "hSphericity", "Sphericity; Sphericity; Entries"},
    {kPoszVsMult, kTH2F, "hPoszVsMult", "Vertex Z vs Multiplicity; V_{Z} (cm); Multiplicity"},
    {kPoszVsCent, kTH2F, "hPoszVsCent", "Vertex Z vs Centrality; V_{Z} (cm); Centrality (%)"},
    {kCentVsMult, kTH2F, "kCentVsMult", "Centrality vs Multiplicity; Centrality (%); Multiplicity"},
    {kMultVsSphericity, kTH2F, "hMultVsSphericity", "Multiplicity vs Sphericity; Multiplicity; Sphericity"},
    {kCentVsSphericity, kTH2F, "kCentVsSphericity", "Centrality vs Sphericity; Centrality (%); Sphericity"},
  }};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
class CollisionHistManager
{
 public:
  /// Destructor
  virtual ~CollisionHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  template <modes::Mode mode>
  void init(HistogramRegistry* registry, std::map<ColHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;
    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(ColAnalysisDir);
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosz, HistTable), GetHistDesc(kPosz, HistTable), GetHistType(kPosz, HistTable), {Specs[kPosz]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMult, HistTable), GetHistDesc(kMult, HistTable), GetHistType(kMult, HistTable), {Specs[kMult]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kCent, HistTable), GetHistDesc(kCent, HistTable), GetHistType(kCent, HistTable), {Specs[kCent]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kSphericity, HistTable), GetHistDesc(kSphericity, HistTable), GetHistType(kSphericity, HistTable), {Specs[kSphericity]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMagField, HistTable), GetHistDesc(kMagField, HistTable), GetHistType(kMagField, HistTable), {Specs[kMagField]});
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {

      std::string qaDir = std::string(ColQaDir);
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPoszVsMult, HistTable), GetHistDesc(kPoszVsMult, HistTable), GetHistType(kPoszVsMult, HistTable), {Specs[kPoszVsMult]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPoszVsCent, HistTable), GetHistDesc(kPoszVsCent, HistTable), GetHistType(kPoszVsCent, HistTable), {Specs[kPoszVsCent]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kCentVsMult, HistTable), GetHistDesc(kCentVsMult, HistTable), GetHistType(kCentVsMult, HistTable), {Specs[kCentVsMult]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kMultVsSphericity, HistTable), GetHistDesc(kMultVsSphericity, HistTable), GetHistType(kMultVsSphericity, HistTable), {Specs[kMultVsSphericity]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kCentVsSphericity, HistTable), GetHistDesc(kCentVsSphericity, HistTable), GetHistType(kCentVsSphericity, HistTable), {Specs[kCentVsSphericity]});
    }
  } // namespace o2::analysis::femtounited

  template <modes::Mode mode, typename T>
  void fill(T const& col)
  {
    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kPosz, HistTable)), col.posZ());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kMult, HistTable)), col.mult());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kCent, HistTable)), col.cent());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kSphericity, HistTable)), col.sphericity());
      mHistogramRegistry->fill(HIST(ColAnalysisDir) + HIST(GetHistName(kMagField, HistTable)), col.magField());
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(GetHistName(kPoszVsMult, HistTable)), col.posZ(), col.mult());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(GetHistName(kPoszVsCent, HistTable)), col.posZ(), col.cent());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(GetHistName(kCentVsMult, HistTable)), col.cent(), col.mult());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(GetHistName(kMultVsSphericity, HistTable)), col.mult(), col.sphericity());
      mHistogramRegistry->fill(HIST(ColQaDir) + HIST(GetHistName(kCentVsSphericity, HistTable)), col.cent(), col.sphericity());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;
}; // namespace femtounitedcolhistmanager
}; // namespace colhistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_COLLISIONHISTMANAGER_H_
