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

/// \file TrackHistmanager.h
/// \brief Femtounited collision histograms

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKHISTMANAGER_H_

#include <array>
#include <vector>
#include <string>
#include <map>

#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoUnited/Core/HistManager.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

using namespace ::o2::framework;

namespace o2::analysis::femtounited
{
namespace trackhistmanager
{
// enum for track histograms
enum TrackHist {
  // kinemtics
  kPt,
  kEta,
  kPhi,
  // qa variables
  kSign,
  kItsCluster,
  kItsClusterIb,
  kTpcCluster,
  kTpcClusterShared,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsItsCluster,
  kPtVsTpcCluster,
  kPtVsTpcClusterShared,
  kTpcClusterVsTpcClusterShared,
  // tpc pid
  kTpcSignal,
  kTpcElectron,
  kTpcPion,
  kTpcKaon,
  kTpcProton,
  kTpcDeuteron,
  kTpcTriton,
  kTpcHelium,
  // tof pid
  kTofBeta,
  kTofElectron,
  kTofPion,
  kTofKaon,
  kTofProton,
  kTofDeuteron,
  kTofTriton,
  kTofHelium,
  // tpc+tof pid
  kTpctofElectron,
  kTpctofPion,
  kTpctofKaon,
  kTpctofProton,
  kTpctofDeuteron,
  kTpctofTriton,
  kTpctofHelium,
  kTrackHistogramLast
};

constexpr std::string_view AnalysisDir = "TrackHistograms/Analysis/";
constexpr std::string_view QaDir = "TrackHistograms/QA/";
constexpr std::string_view PidDir = "TrackHistograms/PID/";

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<Histmanager::HistInfo<TrackHist>, kTrackHistogramLast> HistTable = {
  {{kPt, kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kSign, kTH1F, "hSign", "Sign of charge ; Sign; Entries"},
   {kItsCluster, kTH1F, "hItsCluster", "ITS cluster; ITS cluster; Entries"},
   {kItsClusterIb, kTH1F, "hItsClusterIb", "ITS cluster in inner barrel; ITS IB cluster; Entries"},
   {kTpcCluster, kTH1F, "hTpcCluster", "TPC cluster found; TPC cluster found; Entries"},
   {kTpcClusterShared, kTH1F, "hTpcClusterShared", "TPC cluster shared; TPC cluster shared ; Entries"},
   {kPtVsEta, kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
   {kPtVsPhi, kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
   {kPhiVsEta, kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
   {kPtVsItsCluster, kTH2F, "hPtVsItsCluster", "p_{T} vs ITS cluster; p_{T} (GeV/#it{c}) ; ITS cluster"},
   {kPtVsTpcCluster, kTH2F, "hPtVsTpcCluster", "p_{T} vs TPC cluster found; p_{T} (GeV/#it{c}) ; TPC cluster found"},
   {kPtVsTpcClusterShared, kTH2F, "hPtVsTpcClusterShared", "p_{T} vs TPC cluster shared; p_{T} (GeV/#it{c}) ; TPC cluster shared"},
   {kTpcClusterVsTpcClusterShared, kTH2F, "hTpcClusterVsTpcClusterShared", "TPC cluster found vs TPC cluster shared; TPC cluster found; TPC cluster shared"},
   {kTpcSignal, kTH2F, "hTpcSignal", "TPC Signal; p (GeV/#it{c}) ; TPC Signal"},
   {kTpcElectron, kTH2F, "hTpcPidElectron", "TPC PID Electron; p (GeV/#it{c}) ; n#sigma_{TPC,el}"},
   {kTpcPion, kTH2F, "hTpcPidPion", "TPC PID Pion; p (GeV/#it{c}) ; n#sigma_{TPC,pi}"},
   {kTpcKaon, kTH2F, "hTpcPidKaon", "TPC PID Kaon; p (GeV/#it{c}) ; n#sigma_{TPC,ka}"},
   {kTpcProton, kTH2F, "hTpcPidProton", "TPC PID Proton; p (GeV/#it{c}) ; n#sigma_{TPC,pr}"},
   {kTpcDeuteron, kTH2F, "hTpcPidDeuteron", "TPC PID Deuteron; p (GeV/#it{c}) ; n#sigma_{TPC,de}"},
   {kTpcTriton, kTH2F, "hTpcPidTriton", "TPC PID Triton; p (GeV/#it{c}) ; n#sigma_{TPC,tr}"},
   {kTpcHelium, kTH2F, "hTpcPidHelium", "TPC PID Helium; p (GeV/#it{c}) ; n#sigma_{TPC,he}"},
   {kTofBeta, kTH2F, "hTofBeta", "TOF #beta; p (GeV/#it{c}) ; TOF #beta"},
   {kTofElectron, kTH2F, "hTofPidElectron", "TOF PID Electron; p (GeV/#it{c}) ; n#sigma_{TOF,el}"},
   {kTofPion, kTH2F, "hTofPidPion", "TOF PID Pion; p (GeV/#it{c}) ; n#sigma_{TOF,pi}"},
   {kTofKaon, kTH2F, "hTofPidKaon", "TOF PID Kaon; p (GeV/#it{c}) ; n#sigma_{TOF,ka}"},
   {kTofProton, kTH2F, "hTofPidProton", "TOF PID Proton; p (GeV/#it{c}) ; n#sigma_{TOF,pr}"},
   {kTofDeuteron, kTH2F, "hTofPidDeuteron", "TOF PID Deuteron; p (GeV/#it{c}) ; n#sigma_{TOF,de}"},
   {kTofTriton, kTH2F, "hTofPidTriton", "TOF PID Triton; p (GeV/#it{c}) ; n#sigma_{TOF,tr}"},
   {kTofHelium, kTH2F, "hTofPidHelium", "TOF PID Helium; p (GeV/#it{c}) ; n#sigma_{TOF,he}"},
   {kTpctofElectron, kTH2F, "hTpctofPidElectron", "TOF PID Electron; p (GeV/#it{c}) ; n#sigma_{TOF,el}"},
   {kTpctofPion, kTH2F, "hTpctofPidPion", "TPC+TOF PID Pion; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,pi}^{2}+n#sigma_{TOF,pi}^{2}}"},
   {kTpctofKaon, kTH2F, "hTpctofPidKaon", "TPC+TOF PID Kaon; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,ka}^{2}+n#sigma_{TOF,ka}^{2}}"},
   {kTpctofProton, kTH2F, "hTpctofPidProton", "TPC+TOF PID Proton; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,pr}^{2}+n#sigma_{TOF,pr}^{2}}"},
   {kTpctofDeuteron, kTH2F, "hTpctofPidDeuteron", "TPC+TOF PID Deuteron; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,de}^{2}+n#sigma_{TOF,de}^{2}}"},
   {kTpctofTriton, kTH2F, "hTpctofPidTriton", "TPC+TOF PID Triton; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,tr}^{2}+n#sigma_{TOF,tr}^{2}}"},
   {kTpctofHelium, kTH2F, "hTpctofPidHelium", "TPC+TOF PID Helium; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,he}^{2}+n#sigma_{TOF,he}^{2}}"}}};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
class TrackHistManager
{
 public:
  /// Destructor
  virtual ~TrackHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(HistogramRegistry* registry, std::map<TrackHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt, HistTable), GetHistDesc(kPt, HistTable), GetHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kEta, HistTable), GetHistDesc(kEta, HistTable), GetHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPhi, HistTable), GetHistDesc(kPhi, HistTable), GetHistType(kPhi, HistTable), {Specs[kPhi]});
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {
      std::string qaDir = std::string(QaDir);

      mHistogramRegistry->add(qaDir + GetHistNamev2(kSign, HistTable), GetHistDesc(kSign, HistTable), GetHistType(kSign, HistTable), {Specs[kSign]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kItsCluster, HistTable), GetHistDesc(kItsCluster, HistTable), GetHistType(kItsCluster, HistTable), {Specs[kItsCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kItsClusterIb, HistTable), GetHistDesc(kItsClusterIb, HistTable), GetHistType(kItsClusterIb, HistTable), {Specs[kItsClusterIb]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTpcCluster, HistTable), GetHistDesc(kTpcCluster, HistTable), GetHistType(kTpcCluster, HistTable), {Specs[kTpcCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTpcClusterShared, HistTable), GetHistDesc(kTpcClusterShared, HistTable), GetHistType(kTpcClusterShared, HistTable), {Specs[kTpcClusterShared]});

      // qa 2d
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsEta, HistTable), GetHistDesc(kPtVsEta, HistTable), GetHistType(kPtVsEta, HistTable), {Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsPhi, HistTable), GetHistDesc(kPtVsPhi, HistTable), GetHistType(kPtVsPhi, HistTable), {Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPhiVsEta, HistTable), GetHistDesc(kPhiVsEta, HistTable), GetHistType(kPhiVsEta, HistTable), {Specs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsItsCluster, HistTable), GetHistDesc(kPtVsItsCluster, HistTable), GetHistType(kPtVsItsCluster, HistTable), {Specs[kPtVsItsCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsTpcCluster, HistTable), GetHistDesc(kPtVsTpcCluster, HistTable), GetHistType(kPtVsTpcCluster, HistTable), {Specs[kPtVsTpcCluster]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kPtVsTpcClusterShared, HistTable), GetHistDesc(kPtVsTpcClusterShared, HistTable), GetHistType(kPtVsTpcClusterShared, HistTable), {Specs[kPtVsTpcClusterShared]});
      mHistogramRegistry->add(qaDir + GetHistNamev2(kTpcClusterVsTpcClusterShared, HistTable), GetHistDesc(kTpcClusterVsTpcClusterShared, HistTable), GetHistType(kTpcClusterVsTpcClusterShared, HistTable), {Specs[kTpcClusterVsTpcClusterShared]});

      std::string pidDir = std::string(PidDir);

      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcSignal, HistTable), GetHistDesc(kTpcSignal, HistTable), GetHistType(kTpcSignal, HistTable), {Specs[kTpcSignal]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcElectron, HistTable), GetHistDesc(kTpcElectron, HistTable), GetHistType(kTpcElectron, HistTable), {Specs[kTpcElectron]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcPion, HistTable), GetHistDesc(kTpcPion, HistTable), GetHistType(kTpcPion, HistTable), {Specs[kTpcPion]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcKaon, HistTable), GetHistDesc(kTpcKaon, HistTable), GetHistType(kTpcKaon, HistTable), {Specs[kTpcKaon]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcProton, HistTable), GetHistDesc(kTpcProton, HistTable), GetHistType(kTpcProton, HistTable), {Specs[kTpcProton]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcDeuteron, HistTable), GetHistDesc(kTpcDeuteron, HistTable), GetHistType(kTpcDeuteron, HistTable), {Specs[kTpcDeuteron]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcTriton, HistTable), GetHistDesc(kTpcTriton, HistTable), GetHistType(kTpcTriton, HistTable), {Specs[kTpcTriton]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpcHelium, HistTable), GetHistDesc(kTpcHelium, HistTable), GetHistType(kTpcHelium, HistTable), {Specs[kTpcHelium]});

      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofBeta, HistTable), GetHistDesc(kTofBeta, HistTable), GetHistType(kTofBeta, HistTable), {Specs[kTofBeta]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofElectron, HistTable), GetHistDesc(kTofElectron, HistTable), GetHistType(kTofElectron, HistTable), {Specs[kTofElectron]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofPion, HistTable), GetHistDesc(kTofPion, HistTable), GetHistType(kTofPion, HistTable), {Specs[kTofPion]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofKaon, HistTable), GetHistDesc(kTofKaon, HistTable), GetHistType(kTofKaon, HistTable), {Specs[kTofKaon]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofProton, HistTable), GetHistDesc(kTofProton, HistTable), GetHistType(kTofProton, HistTable), {Specs[kTofProton]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofDeuteron, HistTable), GetHistDesc(kTofDeuteron, HistTable), GetHistType(kTofDeuteron, HistTable), {Specs[kTofDeuteron]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofTriton, HistTable), GetHistDesc(kTofTriton, HistTable), GetHistType(kTofTriton, HistTable), {Specs[kTofTriton]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTofHelium, HistTable), GetHistDesc(kTofHelium, HistTable), GetHistType(kTofHelium, HistTable), {Specs[kTofHelium]});

      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpctofElectron, HistTable), GetHistDesc(kTpctofElectron, HistTable), GetHistType(kTpctofElectron, HistTable), {Specs[kTpctofElectron]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpctofPion, HistTable), GetHistDesc(kTpctofPion, HistTable), GetHistType(kTpctofPion, HistTable), {Specs[kTpctofPion]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpctofKaon, HistTable), GetHistDesc(kTpctofKaon, HistTable), GetHistType(kTpctofKaon, HistTable), {Specs[kTpctofKaon]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpctofProton, HistTable), GetHistDesc(kTpctofProton, HistTable), GetHistType(kTpctofProton, HistTable), {Specs[kTpctofProton]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpctofDeuteron, HistTable), GetHistDesc(kTpctofDeuteron, HistTable), GetHistType(kTpctofDeuteron, HistTable), {Specs[kTpctofDeuteron]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpctofTriton, HistTable), GetHistDesc(kTpctofTriton, HistTable), GetHistType(kTpctofTriton, HistTable), {Specs[kTpctofTriton]});
      mHistogramRegistry->add(pidDir + GetHistNamev2(kTpctofHelium, HistTable), GetHistDesc(kTpctofHelium, HistTable), GetHistType(kTpctofHelium, HistTable), {Specs[kTpctofHelium]});
    }
  }

  template <modes::Mode mode, typename T>
  void fill(T const& track)
  {
    if constexpr (isModeSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPt, HistTable)), track.pt());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kEta, HistTable)), track.eta());
      mHistogramRegistry->fill(HIST(AnalysisDir) + HIST(GetHistName(kPhi, HistTable)), track.phi());
    }

    if constexpr (isModeSet(mode, modes::Mode::kQA)) {
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kSign, HistTable)), static_cast<float>(track.sign()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kItsCluster, HistTable)), static_cast<float>(track.itsNCls()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kItsClusterIb, HistTable)), static_cast<float>(track.itsNClsInnerBarrel()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kTpcCluster, HistTable)), static_cast<float>(track.tpcNClsFound()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kTpcClusterShared, HistTable)), static_cast<float>(track.tpcNClsShared()));

      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsEta, HistTable)), track.pt(), track.eta());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsPhi, HistTable)), track.pt(), track.phi());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPhiVsEta, HistTable)), track.phi(), track.eta());
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsItsCluster, HistTable)), track.pt(), static_cast<float>(track.itsNCls()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsTpcCluster, HistTable)), track.pt(), static_cast<float>(track.tpcNClsFound()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kPtVsTpcClusterShared, HistTable)), track.pt(), static_cast<float>(track.tpcNClsShared()));
      mHistogramRegistry->fill(HIST(QaDir) + HIST(GetHistName(kTpcClusterVsTpcClusterShared, HistTable)), static_cast<float>(track.tpcNClsFound()), static_cast<float>(track.tpcNClsShared()));

      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcSignal, HistTable)), track.p(), track.tpcSignal());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcElectron, HistTable)), track.p(), track.tpcNSigmaEl());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcPion, HistTable)), track.p(), track.tpcNSigmaPi());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcKaon, HistTable)), track.p(), track.tpcNSigmaKa());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcProton, HistTable)), track.p(), track.tpcNSigmaPr());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcDeuteron, HistTable)), track.p(), track.tpcNSigmaDe());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcTriton, HistTable)), track.p(), track.tpcNSigmaTr());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpcHelium, HistTable)), track.p(), track.tpcNSigmaHe());

      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofBeta, HistTable)), track.p(), track.tofBeta());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofElectron, HistTable)), track.p(), track.tofNSigmaEl());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofPion, HistTable)), track.p(), track.tofNSigmaPi());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofKaon, HistTable)), track.p(), track.tofNSigmaKa());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofProton, HistTable)), track.p(), track.tofNSigmaPr());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofDeuteron, HistTable)), track.p(), track.tofNSigmaDe());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofTriton, HistTable)), track.p(), track.tofNSigmaTr());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTofHelium, HistTable)), track.p(), track.tofNSigmaHe());

      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpctofElectron, HistTable)), track.p(), track.tpctofNSigmaEl());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpctofPion, HistTable)), track.p(), track.tpctofNSigmaPi());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpctofKaon, HistTable)), track.p(), track.tpctofNSigmaKa());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpctofProton, HistTable)), track.p(), track.tpctofNSigmaPr());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpctofDeuteron, HistTable)), track.p(), track.tpctofNSigmaDe());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpctofTriton, HistTable)), track.p(), track.tpctofNSigmaTr());
      mHistogramRegistry->fill(HIST(PidDir) + HIST(GetHistName(kTpctofHelium, HistTable)), track.p(), track.tpctofNSigmaHe());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;
}; // namespace femtounitedcolhistmanager
}; // namespace trackhistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_TRACKHISTMANAGER_H_
