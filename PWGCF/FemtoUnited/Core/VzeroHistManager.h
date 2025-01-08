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

/// \file VzeroHistManager.h
/// \brief histogram manager for vzero histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_VZEROHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_VZEROHISTMANAGER_H_

#include <array>
#include <vector>
#include <string>
#include <map>

#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoUnited/Core/HistManager.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"
#include "PWGCF/FemtoUnited/Utils/FemtoUtils.h"

namespace o2::analysis::femtounited
{
namespace vzerohistmanager
{
// enum for track histograms
enum VzeroHist {
  // analysis
  kPt,
  kEta,
  kPhi,
  kMass,
  kPosDauPt,
  kPosDauEta,
  kPosDauPhi,
  kNegDauPt,
  kNegDauEta,
  kNegDauPhi,
  // qa variables
  kSign,
  kDecayDauDca,
  kDecayVtxX,
  kDecayVtxY,
  kDecayVtxZ,
  kTransRadius,
  kKaonMass,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  // qa for daughters
  kPosDaughTpcCluster,
  kPosDaughPtVsDcaxy,
  kPosDaughPtVsDcaz,
  kPosDaughPtVsDca,
  kPosDaughTpcNsigma,
  kNegDaughTpcCluster,
  kNegDaughPtVsDcaxy,
  kNegDaughPtVsDcaz,
  kNegDaughPtVsDca,
  kNegDaughTpcNsigma,
  kVeroHistLast
};

constexpr std::string_view AnalysisDir = "VzeroHistograms/Analysis/";
constexpr std::string_view QaDir = "VzeroHistograms/QA/";
constexpr std::string_view PidDir = "VzeroHistograms/Daughters/";

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<VzeroHist>, kVeroHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kMass, o2::framework::kTH1F, "hMass", "Invariant mass; m_{Inv} (GeV/#it{c}^{2}; Entries"},
   {kPosDauPt, o2::framework::kTH1F, "hPosDauPt", "Transverse Momentum of positive daughter; p_{T} (GeV/#it{c}); Entries"},
   {kPosDauEta, o2::framework::kTH1F, "hPosDauEta", "Pseudorapdity of positive daughter; #eta; Entries"},
   {kPosDauPhi, o2::framework::kTH1F, "hPosDauPhi", "Azimuthal angle of positive daughter; #varphi; Entries"},
   {kNegDauPt, o2::framework::kTH1F, "hNegDauPt", "Transverse Momentum of negative daughter; p_{T} (GeV/#it{c}); Entries"},
   {kNegDauEta, o2::framework::kTH1F, "hNegDauEta", "Pseudorapdity of negative daughter; #eta; Entries"},
   {kNegDauPhi, o2::framework::kTH1F, "hNegDauPhi", "Azimuthal angle of negative daughter; #varphi; Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign (+1 particle/-1 antiparticle) ; Sign; Entries"},
   {kDecayDauDca, o2::framework::kTH1F, "hDauDca", "Daughter DCA at decay vertex ; DCA_{Decay vertex} (cm); Entries"},
   {kDecayVtxX, o2::framework::kTH1F, "hDecayVtxX", "X coordinate of decay vertex ; DV_{X} (cm); Entries"},
   {kDecayVtxY, o2::framework::kTH1F, "hDecayVtxY", "Y coordinate of decay vertex ; DV_{Y} (cm); Entries"},
   {kDecayVtxZ, o2::framework::kTH1F, "hDecayVtxZ", "Z coordinate of decay vertex ; DV_{Z} (cm); Entries"},
   {kTransRadius, o2::framework::kTH1F, "hTransRadius", "Tranverse radius ; r_{xy} (cm); Entries"},
   {kKaonMass, o2::framework::kTH1F, "hKaonMass", "Kaon mass hypothesis ; m_{K} (GeV/#it{c}^{2}); Entries"},
   {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPosDaughTpcCluster, o2::framework::kTH1F, "hPosDaughTpcCluster", "TPC cluster found (daughter^{+}) ; TPC cluster found, Entries"},
   {kPosDaughPtVsDcaxy, o2::framework::kTH2F, "hPosDauPtVsDcaxy", "p_{T} vs DCA_{XY} (daughter +) ; p_{T} (GeV/#it{c}) ; DCA_{XY} (cm)"},
   {kPosDaughPtVsDcaz, o2::framework::kTH2F, "hPosDauPtVsDcaz", "p_{T} vs DCA_{Z} (daughter +) ; p_{T} (GeV/#it{c}) ; DCA_{Z}"},
   {kPosDaughPtVsDca, o2::framework::kTH2F, "hPosDauPtVsDca", "p_{T} vs DCA (daughter +) ; p_{T} (GeV/#it{c}) ; DCA"},
   {kPosDaughTpcNsigma, o2::framework::kTH2F, "hPosDauTpcNsigma", "TPC PID (daughter +) ; p (GeV/#it{c}) ; n#sigma_{TPC}"},
   {kNegDaughTpcCluster, o2::framework::kTH1F, "hNegDAugTpcCluster", "TPC cluster found (daughter -) ; TPC cluster found; Entries"},
   {kNegDaughPtVsDcaxy, o2::framework::kTH2F, "hNegDauPtVsDcaxy", "p_{T} vs DCA_{XY} (daughter -) ; p_{T} (GeV/#it{c}) ; DCA_{XY} (cm)"},
   {kNegDaughPtVsDcaz, o2::framework::kTH2F, "hNegDauPtVsDcaz", "p_{T} vs DCA_{Z} (daughter -) ; p_{T} (GeV/#it{c}) ; DCA_{Z}"},
   {kNegDaughPtVsDca, o2::framework::kTH2F, "hNegDauPtVsDca", "p_{T} vs DCA (daughter -) ; p_{T} (GeV/#it{c}) ; DCA"},
   {kNegDaughTpcNsigma, o2::framework::kTH2F, "hNegDaughTpcNsigma", "TPC PID (daughter -) ; p (GeV/#it{c}) ; n#sigma_{TPC}"}}};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
class VzeroHistManager
{
 public:
  /// Destructor
  virtual ~VzeroHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<VzeroHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMass, HistTable), GetHistDesc(kMass, HistTable), GetHistType(kMass, HistTable), {Specs[kMass]});

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosDauPt, HistTable), GetHistDesc(kPosDauPt, HistTable), GetHistType(kPosDauPt, HistTable), {Specs[kPosDauPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosDauEta, HistTable), GetHistDesc(kPosDauEta, HistTable), GetHistType(kPosDauEta, HistTable), {Specs[kPosDauEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPosDauPhi, HistTable), GetHistDesc(kPosDauPhi, HistTable), GetHistType(kPosDauPhi, HistTable), {Specs[kPosDauPhi]});

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kNegDauPt, HistTable), GetHistDesc(kNegDauPt, HistTable), GetHistType(kNegDauPt, HistTable), {Specs[kNegDauPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kNegDauEta, HistTable), GetHistDesc(kNegDauEta, HistTable), GetHistType(kNegDauEta, HistTable), {Specs[kNegDauEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kNegDauPhi, HistTable), GetHistDesc(kNegDauPhi, HistTable), GetHistType(kNegDauPhi, HistTable), {Specs[kNegDauPhi]});
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {
      std::string qaDir = std::string(QaDir);

      mHistogramRegistry->add(qaDir + GetHistNamev2(kSign, HistTable), GetHistDesc(kSign, HistTable), GetHistType(kSign, HistTable), {Specs[kSign]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayDauDca, HistTable), GetHistDesc(kDecayDauDca, HistTable), GetHistType(kDecayDauDca, HistTable), {Specs[kDecayDauDca]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxX, HistTable), GetHistDesc(kDecayVtxX, HistTable), GetHistType(kDecayVtxX, HistTable), {Specs[kDecayVtxX]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxY, HistTable), GetHistDesc(kDecayVtxY, HistTable), GetHistType(kDecayVtxY, HistTable), {Specs[kDecayVtxY]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kDecayVtxZ, HistTable), GetHistDesc(kDecayVtxZ, HistTable), GetHistType(kDecayVtxZ, HistTable), {Specs[kDecayVtxZ]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTransRadius, HistTable), GetHistDesc(kTransRadius, HistTable), GetHistType(kTransRadius, HistTable), {Specs[kTransRadius]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kKaonMass, HistTable), GetHistDesc(kKaonMass, HistTable), GetHistType(kKaonMass, HistTable), {Specs[kKaonMass]});

      // qa 2d
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, HistTable), GetHistDesc(kPtVsEta, HistTable), GetHistType(kPtVsEta, HistTable), {Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, HistTable), GetHistDesc(kPtVsPhi, HistTable), GetHistType(kPtVsPhi, HistTable), {Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, HistTable), GetHistDesc(kPhiVsEta, HistTable), GetHistType(kPhiVsEta, HistTable), {Specs[kPhiVsEta]});

      // qa daughters
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDaughTpcCluster, HistTable), GetHistDesc(kPosDaughTpcCluster, HistTable), GetHistType(kPosDaughTpcCluster, HistTable), {Specs[kPosDaughTpcCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDaughTpcNsigma, HistTable), GetHistDesc(kPosDaughTpcNsigma, HistTable), GetHistType(kPosDaughTpcNsigma, HistTable), {Specs[kPosDaughTpcNsigma]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDaughPtVsDcaxy, HistTable), GetHistDesc(kPosDaughPtVsDcaxy, HistTable), GetHistType(kPosDaughPtVsDcaxy, HistTable), {Specs[kPosDaughPtVsDcaxy]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDaughPtVsDcaz, HistTable), GetHistDesc(kPosDaughPtVsDcaz, HistTable), GetHistType(kPosDaughPtVsDcaz, HistTable), {Specs[kPosDaughPtVsDcaz]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPosDaughPtVsDca, HistTable), GetHistDesc(kPosDaughPtVsDca, HistTable), GetHistType(kPosDaughPtVsDca, HistTable), {Specs[kPosDaughPtVsDca]});

      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDaughTpcCluster, HistTable), GetHistDesc(kNegDaughTpcCluster, HistTable), GetHistType(kNegDaughTpcCluster, HistTable), {Specs[kNegDaughTpcCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDaughPtVsDcaxy, HistTable), GetHistDesc(kNegDaughPtVsDcaxy, HistTable), GetHistType(kNegDaughPtVsDcaxy, HistTable), {Specs[kNegDaughPtVsDcaxy]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDaughPtVsDcaz, HistTable), GetHistDesc(kNegDaughPtVsDcaz, HistTable), GetHistType(kNegDaughPtVsDcaz, HistTable), {Specs[kNegDaughPtVsDcaz]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDaughPtVsDca, HistTable), GetHistDesc(kNegDaughPtVsDca, HistTable), GetHistType(kNegDaughPtVsDca, HistTable), {Specs[kNegDaughPtVsDca]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kNegDaughTpcNsigma, HistTable), GetHistDesc(kNegDaughTpcNsigma, HistTable), GetHistType(kNegDaughTpcNsigma, HistTable), {Specs[kNegDaughTpcNsigma]});
    }
  }

  template <modes::Mode mode, typename T>
  void fill(T const& v0)
  {
    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), v0.pt());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), v0.eta());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), v0.phi());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kMass, HistTable)), v0.vzeroMass());

      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPosDauPt, HistTable)), v0.posDauPt());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPosDauEta, HistTable)), v0.posDauEta());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPosDauPhi, HistTable)), v0.posDauPhi());

      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kNegDauPt, HistTable)), v0.negDauPt());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kNegDauEta, HistTable)), v0.negDauEta());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kNegDauPhi, HistTable)), v0.negDauPhi());
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kSign, HistTable)), static_cast<float>(v0.sign()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayDauDca, HistTable)), v0.dauDCA());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayVtxX, HistTable)), v0.decayVtxX());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayVtxY, HistTable)), v0.decayVtxY());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kDecayVtxZ, HistTable)), v0.decayVtxZ());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kTransRadius, HistTable)), v0.transRadius());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kKaonMass, HistTable)), v0.kaonMass());

      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), v0.pt(), v0.eta());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), v0.pt(), v0.phi());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), v0.phi(), v0.eta());

      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDaughTpcCluster, HistTable)), v0.posDauTpcNclsFound());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDaughPtVsDcaxy, HistTable)), v0.posDauPt(), v0.posDauDcaxy());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDaughPtVsDcaz, HistTable)), v0.posDauPt(), v0.posDauDcaz());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDaughPtVsDca, HistTable)), v0.posDauP(), utils::geometricMean(v0.posDauDcaxy(), v0.posDauDcaz()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPosDaughTpcNsigma, HistTable)), v0.posDauP(), v0.posDauTpcNsigma());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDaughTpcCluster, HistTable)), v0.negDauTpcNclsFound());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDaughPtVsDcaxy, HistTable)), v0.negDauPt(), v0.negDauDcaxy());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDaughPtVsDcaz, HistTable)), v0.negDauPt(), v0.negDauDcaz());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDaughPtVsDca, HistTable)), v0.negDauP(), utils::geometricMean(v0.negDauDcaxy(), v0.negDauDcaz()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kNegDaughTpcNsigma, HistTable)), v0.negDauP(), v0.negDauTpcNsigma());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
};
}; // namespace vzerohistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_VZEROHISTMANAGER_H_
