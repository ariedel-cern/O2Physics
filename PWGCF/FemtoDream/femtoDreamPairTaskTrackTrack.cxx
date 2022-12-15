// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamPairTaskTrackTrack.cxx
/// \brief Tasks that reads the track tables used for the pairing and builds pairs of two tracks
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "FemtoDreamParticleHisto.h"
#include "FemtoDreamEventHisto.h"
#include "FemtoDreamPairCleaner.h"
#include "FemtoDreamContainer.h"
#include "FemtoDreamDetaDphiStar.h"
#include "FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
static constexpr int nPart = 2;
static constexpr int nCuts = 5;
static const std::vector<std::string> partNames{"PartOne", "PartTwo"};
static const std::vector<std::string> cutNames{"MaxPt", "PIDthr", "nSigmaTPC", "nSigmaTPCTOF", "MaxP"};
static const float cutsTable[nPart][nCuts]{
  {4.05f, 1.f, 3.f, 3.f, 100.f},
  {4.05f, 1.f, 3.f, 3.f, 100.f}};
} // namespace

struct femtoDreamPairTaskTrackTrack {

  /// Particle selection part

  /// Table for both particles
  Configurable<LabeledArray<float>> cfgCutTable{"cfgCutTable", {cutsTable[0], nPart, nCuts, partNames, cutNames}, "Particle selections"};
  Configurable<int> cfgNspecies{"ccfgNspecies", 4, "Number of particle spieces with PID info"};
  Configurable<std::vector<float>> ConfPIDnSigmaMax{"ConfPIDnSigmaMax", std::vector<float>{3.5f, 3.f, 2.5f}, "This configurable needs to be the same as the one used in the producer task"};
  Configurable<bool> isMonteCarlo{"isMonteCarlo", false, "Enable additional Histogramms in the case of a MonteCarlo Run"};

  /// Particle 1
  Configurable<int> ConfPDGCodePartOne{"ConfPDGCodePartOne", 2212, "Particle 1 - PDG code"};
  Configurable<uint32_t> ConfCutPartOne{"ConfCutPartOne", 5542474, "Particle 1 - Selection bit from cutCulator"};
  Configurable<std::vector<int>> ConfPIDPartOne{"ConfPIDPartOne", std::vector<int>{2}, "Particle 1 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>int>>
  ConfigurableAxis CfgpTBinsPartOne{"CfgpTBinsPartOne", {20, 0.5, 4.05}, "Particle 1 - binning pT in 2D plots"};
  ConfigurableAxis CfgDCAxyBinsPartOne{"CfgDCAxyBinsPartOne", {500, -0.5, 0.5}, "Particle 1 - binning DCAxy"};

  /// Partition for particle 1
  Partition<aod::FemtoDreamParticles> partsOne = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);
  Partition<Join<aod::FemtoDreamParticles, aod::FemtoDreamParticlesMC>> partsOneMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartOne) == ConfCutPartOne);

  /// Histogramming for particle 1
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 1> trackHistoPartOne;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 3> trackHistoPartOneMC;

  /// Particle 2
  Configurable<bool> ConfIsSame{"ConfIsSame", false, "Pairs of the same particle"};
  Configurable<int> ConfPDGCodePartTwo{"ConfPDGCodePartTwo", 2212, "Particle 2 - PDG code"};
  Configurable<uint32_t> ConfCutPartTwo{"ConfCutPartTwo", 5542474, "Particle 2 - Selection bit"};
  Configurable<std::vector<int>> ConfPIDPartTwo{"ConfPIDPartTwo", std::vector<int>{2}, "Particle 2 - Read from cutCulator"}; // we also need the possibility to specify whether the bit is true/false ->std>>vector<std::pair<int, int>>
  ConfigurableAxis CfgpTBinsPartTwo{"CfgpTBinsPartTwo", {20, 0.5, 4.05}, "Particle 2 - binning pT in 2D plots"};
  ConfigurableAxis CfgDCAxyBinsPartTwo{"CfgDCAxyBinsPartTwo", {500, -0.5, 0.5}, "Particle 2 - binning DCAxy"};

  /// Partition for particle 2
  Partition<aod::FemtoDreamParticles> partsTwo = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);
  Partition<Join<aod::FemtoDreamParticles, aod::FemtoDreamParticlesMC>> partsTwoMC = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) && ((aod::femtodreamparticle::cut & ConfCutPartTwo) == ConfCutPartTwo);

  /// Histogramming for particle 2
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 2> trackHistoPartTwo;
  FemtoDreamParticleHisto<aod::femtodreamparticle::ParticleType::kTrack, 4> trackHistoPartTwoMC;

  /// Histogramming for Event
  FemtoDreamEventHisto eventHisto;

  /// The configurables need to be passed to an std::vector
  std::vector<int> vPIDPartOne, vPIDPartTwo;
  std::vector<float> kNsigma;

  /// Correlation part
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"}; // \todo to be obtained from the hash task
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  ColumnBinningPolicy<aod::collision::PosZ, aod::femtodreamcollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

  ConfigurableAxis CfgkstarBins{"CfgkstarBins", {1500, 0., 6.}, "binning kstar"};
  ConfigurableAxis CfgkTBins{"CfgkTBins", {150, 0., 9.}, "binning kT"};
  ConfigurableAxis CfgmTBins{"CfgmTBins", {225, 0., 7.5}, "binning mT"};
  Configurable<int> ConfNEventsMix{"ConfNEventsMix", 5, "Number of events for mixing"};
  Configurable<bool> ConfIsCPR{"ConfIsCPR", true, "Close Pair Rejection"};
  Configurable<bool> ConfCPRPlotPerRadii{"ConfCPRPlotPerRadii", false, "Plot CPR per radii"};

  FemtoDreamContainer<femtoDreamContainer::EventType::same, femtoDreamContainer::Observable::kstar> sameEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::sameMC, femtoDreamContainer::Observable::kstar> sameEventContMC;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixed, femtoDreamContainer::Observable::kstar> mixedEventCont;
  FemtoDreamContainer<femtoDreamContainer::EventType::mixedMC, femtoDreamContainer::Observable::kstar> mixedEventContMC;

  FemtoDreamPairCleaner<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCleaner;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> pairCloseRejection;
  /// Histogram output
  HistogramRegistry qaRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry MixQaRegistry{"MixQaRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&qaRegistry);
    trackHistoPartOne.init(&qaRegistry, CfgpTBinsPartOne, CfgDCAxyBinsPartOne);
    if (isMonteCarlo) {
      trackHistoPartOneMC.initMC(&qaRegistry, CfgpTBinsPartOne, CfgDCAxyBinsPartOne);
    }
    if (!ConfIsSame) {
      trackHistoPartTwo.init(&qaRegistry, CfgpTBinsPartTwo, CfgDCAxyBinsPartTwo);
      if (isMonteCarlo) {
        trackHistoPartTwoMC.initMC(&qaRegistry, CfgpTBinsPartTwo, CfgDCAxyBinsPartTwo);
      }
    }

    MixQaRegistry.add("MixingQA/hSECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});
    MixQaRegistry.add("MixingQA/hMECollisionBins", ";bin;Entries", kTH1F, {{120, -0.5, 119.5}});

    sameEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
    sameEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);

    if (isMonteCarlo) {
      sameEventContMC.initMC(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
      sameEventContMC.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    }

    mixedEventCont.init(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
    mixedEventCont.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);

    if (isMonteCarlo) {
      mixedEventContMC.initMC(&resultRegistry, CfgkstarBins, CfgMultBins, CfgkTBins, CfgmTBins);
      mixedEventContMC.setPDGCodes(ConfPDGCodePartOne, ConfPDGCodePartTwo);
    }

    pairCleaner.init(&qaRegistry);
    if (ConfIsCPR) {
      pairCloseRejection.init(&resultRegistry, &qaRegistry, 0.01, 0.01, ConfCPRPlotPerRadii); /// \todo add config for Δη and ΔΦ cut values
    }

    vPIDPartOne = ConfPIDPartOne;
    vPIDPartTwo = ConfPIDPartTwo;
    kNsigma = ConfPIDnSigmaMax;
  }

  /// This function processes the same event and takes care of all the histogramming
  /// \todo the trivial loops over the tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  template <bool isMC, typename CollisionType, typename PartType, typename SlicedGroupType>
  void doSameEvent(CollisionType const& col, PartType const& parts, SlicedGroupType const& groupPartsOne, SlicedGroupType const& groupPartsTwo)
  {

    MixQaRegistry.fill(HIST("MixingQA/hSECollisionBins"), colBinning.getBin({col.posZ(), col.multNtrPV()}));

    const auto& magFieldTesla = col.magField();

    const int multCol = col.multNtrPV();
    eventHisto.fillQA(col);
    /// Histogramming same event
    for (auto& part : groupPartsOne) {
      if (part.p() > cfgCutTable->get("PartOne", "MaxP") || part.pt() > cfgCutTable->get("PartOne", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(part.pidcut(),
                             part.p(),
                             cfgCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartOne", "nSigmaTPC"),
                             cfgCutTable->get("PartOne", "nSigmaTPCTOF"))) {
        continue;
      }
      trackHistoPartOne.fillQA(part);
      if constexpr (isMC) {
        trackHistoPartOneMC.fillQAMC(part);
      }
    }
    if (!ConfIsSame) {
      for (auto& part : groupPartsTwo) {
        if (part.p() > cfgCutTable->get("PartTwo", "MaxP") || part.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
          continue;
        }
        if (!isFullPIDSelected(part.pidcut(),
                               part.p(),
                               cfgCutTable->get("PartTwo", "PIDthr"),
                               vPIDPartTwo,
                               cfgNspecies,
                               kNsigma,
                               cfgCutTable->get("PartTwo", "nSigmaTPC"),
                               cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
          continue;
        }
        trackHistoPartTwo.fillQA(part);
        if constexpr (isMC) {
          trackHistoPartTwoMC.fillQAMC(part);
        }
      }
    }
    /// Now build the combinations
    for (auto& [p1, p2] : combinations(groupPartsOne, groupPartsTwo)) {
      if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.p() > cfgCutTable->get("PartTwo", "MaxP") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(),
                             p1.p(),
                             cfgCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartOne", "nSigmaTPC"),
                             cfgCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidcut(),
                             p2.p(),
                             cfgCutTable->get("PartTwo", "PIDthr"),
                             vPIDPartTwo,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartTwo", "nSigmaTPC"),
                             cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }

      // track cleaning
      if (!pairCleaner.isCleanPair(p1, p2, parts)) {
        continue;
      }
      sameEventCont.setPair(p1, p2, multCol);
      if constexpr (isMC) {
        sameEventContMC.setPairMC(p1, p2, multCol);
      }
    }
  }

  void processSameEvent(o2::aod::FemtoDreamCollision& col,
                        o2::aod::FemtoDreamParticles& parts)
  {
    auto slicedPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto slicedPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    doSameEvent<false>(col, parts, slicedPartsOne, slicedPartsTwo);
    // doSameEvent<false>(col, parts);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEvent, "Enable processing same event", true);

  void processSameEventMC(o2::aod::FemtoDreamCollision& col,
                          Join<aod::FemtoDreamParticles, aod::FemtoDreamParticlesMC>& parts)
  {
    auto slicedPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto slicedPartsTwo = partsTwoMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    doSameEvent<true>(col, parts, slicedPartsOne, slicedPartsTwo);
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processSameEventMC, "Enable processing same event for Monte Carlo", false);

  /// This function processes the mixed event
  /// \todo the trivial loops over the collisions and tracks should be factored out since they will be common to all combinations of T-T, T-V0, V0-V0, ...
  template <bool isMC, typename PartType, typename SlicedGroupType>
  void doMixedEvent(float const& magFieldTesla, float const& multPV, PartType const& parts, SlicedGroupType const& groupPartsOne, SlicedGroupType const& groupPartsTwo)
  {

    /// \todo before mixing we should check whether both collisions contain a pair of particles!
    // if (partsOne.size() == 0 || nPart2Evt1 == 0 || nPart1Evt2 == 0 || partsTwo.size() == 0 ) continue;

    for (auto& [p1, p2] : combinations(CombinationsFullIndexPolicy(groupPartsOne, groupPartsTwo))) {
      if (p1.p() > cfgCutTable->get("PartOne", "MaxP") || p1.pt() > cfgCutTable->get("PartOne", "MaxPt") || p2.p() > cfgCutTable->get("PartTwo", "MaxP") || p2.pt() > cfgCutTable->get("PartTwo", "MaxPt")) {
        continue;
      }
      if (!isFullPIDSelected(p1.pidcut(),
                             p1.p(),
                             cfgCutTable->get("PartOne", "PIDthr"),
                             vPIDPartOne,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartOne", "nSigmaTPC"),
                             cfgCutTable->get("PartOne", "nSigmaTPCTOF")) ||
          !isFullPIDSelected(p2.pidcut(),
                             p2.p(),
                             cfgCutTable->get("PartTwo", "PIDthr"),
                             vPIDPartTwo,
                             cfgNspecies,
                             kNsigma,
                             cfgCutTable->get("PartTwo", "nSigmaTPC"),
                             cfgCutTable->get("PartTwo", "nSigmaTPCTOF"))) {
        continue;
      }

      if (ConfIsCPR) {
        if (pairCloseRejection.isClosePair(p1, p2, parts, magFieldTesla)) {
          continue;
        }
      }
      mixedEventCont.setPair(p1, p2, multPV);
      if constexpr (isMC) {
        mixedEventContMC.setPairMC(p1, p2, multPV);
      }
    }
  }

  void processMixedEvent(o2::aod::FemtoDreamCollisions& cols,
                         o2::aod::FemtoDreamParticles& parts)
  {
    for (auto& [col1, col2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({col1.posZ(), col1.multNtrPV()}));

      auto slicedPartsOne = partsOne->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col1.globalIndex());
      auto slicedPartsTwo = partsTwo->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col2.globalIndex());

      const auto& magFieldTesla1 = col1.magField();
      const auto& magFieldTesla2 = col2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      doMixedEvent<false>(magFieldTesla1, col1.multNtrPV(), parts, slicedPartsOne, slicedPartsTwo);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEvent, "Enable processing mixed events", true);

  void processMixedEventMC(o2::aod::FemtoDreamCollisions& cols,
                           Join<aod::FemtoDreamParticles, aod::FemtoDreamParticlesMC>& parts)
  {
    for (auto& [col1, col2] : soa::selfCombinations(colBinning, 5, -1, cols, cols)) {

      MixQaRegistry.fill(HIST("MixingQA/hMECollisionBins"), colBinning.getBin({col1.posZ(), col1.multNtrPV()}));

      auto slicedPartsOne = partsOneMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col1.globalIndex());
      auto slicedPartsTwo = partsTwoMC->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col2.globalIndex());

      const auto& magFieldTesla1 = col1.magField();
      const auto& magFieldTesla2 = col2.magField();

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      doMixedEvent<true>(magFieldTesla1, col1.multNtrPV(), parts, slicedPartsOne, slicedPartsTwo);
    }
  }
  PROCESS_SWITCH(femtoDreamPairTaskTrackTrack, processMixedEventMC, "Enable processing mixed events MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamPairTaskTrackTrack>(cfgc),
  };
  return workflow;
}
