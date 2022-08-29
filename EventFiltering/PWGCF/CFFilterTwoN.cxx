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

/// \file CFFilterTwoN.cxx
/// \brief Selection of events with different kind of pairs for femtoscopic studies
///
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include <Framework/Configurable.h>
#include <TMath.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "../filterTables.h"
#include "../../PWGCF/FemtoDream/FemtoUtils.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/FemtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/FemtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/FemtoDreamContainer.h"
#include "PWGCF/FemtoDream/FemtoDreamMath.h"
#include "PWGCF/FemtoDream/FemtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/FemtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/FemtoDreamContainer.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonConstants/PhysicsConstants.h"

#include <vector>

namespace
{

enum kCFTwoBodyTriggers {
  kPD, //=0
  kLD, //=1
  kLAST_CFTwoBodyTriggers
};

static const std::vector<std::string> CfTriggerNames{"kPD", "kLD"};
static constexpr uint8_t Track = 0;
static constexpr uint8_t V0 = 1; // V0
// static constexpr uint8_t V0Daughter = 2; // V0  daughters

static constexpr uint32_t kSignMinusMask = 1;
static constexpr uint32_t kSignPlusMask = 2;

// static constexpr uint32_t knSigmaProton = 48;
static constexpr uint32_t kValue0 = 0;

} // namespace

namespace o2::aod
{
using FullCollision = soa::Join<aod::Collisions,
                                aod::EvSels,
                                aod::Mults>::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtoDream;

struct CFFilterTwoN {

  Produces<aod::CFFiltersTwoN> tags;

  Configurable<float> confKstarTriggerLimit{"KstarTriggerLimitUpper", 1.0f, "Kstar limit for selection"};
  Configurable<int> KstarTrigger{"KstarTrigger", 0, "Choice which trigger to run"};

  Configurable<float> confProtonPtMin{"ProtonPtMin", 0.5, "Minimal Pt for Protons"};
  Configurable<float> confProtonPtMax{"ProtonPtMax", 4.0, "Maximal Pt for Protons"};
  Configurable<float> confPIDThresholdProton{"PIDThresholdProton", 0.75f, "P threshold for TPC/TPC&TOF selection (Protons only)"};

  Configurable<float> confDeuteronPtMin{"DeuteronPtMin", 0.5, "Minimal Pt for Deuterons"};
  Configurable<float> confDeuteronPtMax{"DeuteronPtMax", 1.4, "Maximal Pt for Deuterons"};

  Configurable<float> confnSigmaAccetance{"PIDAcceptance", 3., "nSigma for accepting Protons and Deuterons"};
  Configurable<float> confPIDRejection{"PIDRejection", 3.5, "nSigma for rejection bogus Deuterons"};

  Configurable<std::vector<float>> confnSigmaReduced{"ReducedNSigma", {1.0, 0.5, 0.3}, "reduce nSigma as a function of pt for Deuteron acceptance"};
  Configurable<std::vector<float>> confnSigmaReducedPt{"ReducedNSigmaPt", {1.1, 1.2, 1.3}, "pt values for reduce nSigma for Deuteron acceptance"};

  Configurable<float> ldeltaPhiMax{"ldeltaPhiMax", 0.010, "Max limit of delta phi"};
  Configurable<float> ldeltaEtaMax{"ldeltaEtaMax", 0.010, "Max limit of delta eta"};

  // obtain particle candidates of protons, deuterons as well as antiprotons and antideuterons
  Partition<o2::aod::FemtoDreamParticles> partPD = (o2::aod::femtodreamparticle::partType == Track) &&
                                                   ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partAntiPD = (o2::aod::femtodreamparticle::partType == Track) &&
                                                       ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);
  // obtain lambdas and antilambdas
  Partition<o2::aod::FemtoDreamParticles> partL = (o2::aod::femtodreamparticle::partType == V0) &&
                                                  ((o2::aod::femtodreamparticle::cut & kSignPlusMask) > kValue0);
  Partition<o2::aod::FemtoDreamParticles> partAntiL = (o2::aod::femtodreamparticle::partType == V0) &&
                                                      ((o2::aod::femtodreamparticle::cut & kSignMinusMask) > kValue0);

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryQA{"registryQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  // containers for close pair rejection
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kTrack> closePairRejectionTT;
  FemtoDreamDetaDphiStar<aod::femtodreamparticle::ParticleType::kTrack, aod::femtodreamparticle::ParticleType::kV0> closePairRejectionTV0;

  bool SelectParticlePID(aod::femtodreamparticle::cutContainerType const& pidCut, int vSpecies, float momentum, float transverseMomentum)
  {
    std::vector<float> vNsigma = {3.5, 3., 2.5, 1., 0.5, 0.3};
    bool pidSelection = false;
    if (vSpecies == o2::track::PID::Proton) {
      // use momentum dependend (TPC or TPC&TOF) pid selection for protons
      pidSelection = isFullPIDSelected(pidCut,
                                       momentum,
                                       confPIDThresholdProton.value,
                                       std::vector<int>{2},
                                       4.,
                                       vNsigma,
                                       confnSigmaAccetance.value,
                                       confnSigmaAccetance);
    } else if (vSpecies == o2::track::PID::Deuteron) {
      // use additional rejection for deuterons
      if (confPIDRejection.value > 0.) {
        // add additional rejections for deuterons if the paritcle could also be a electron, pion or proton
        if (!isPIDSelected(pidCut, std::vector<int>{0}, 4, confPIDRejection.value, vNsigma, kDetector::kTPC) &&
            !isPIDSelected(pidCut, std::vector<int>{1}, 4, confPIDRejection.value, vNsigma, kDetector::kTPC) &&
            !isPIDSelected(pidCut, std::vector<int>{2}, 4, confPIDRejection.value, vNsigma, kDetector::kTPC)) {
          pidSelection = isPIDSelected(pidCut, std::vector<int>{3}, 4, confnSigmaAccetance.value, vNsigma, kDetector::kTPC);
        }
      } else if (!confnSigmaReduced.value.empty()) {
        if (transverseMomentum > confnSigmaReducedPt.value.at(0)) {
          for (std::size_t i = 1; i < confnSigmaReduced.value.size(); i++) {
            if (transverseMomentum < confnSigmaReducedPt.value.at(i)) {
              pidSelection = isPIDSelected(pidCut, std::vector<int>{3}, 4, confnSigmaReduced.value.at(i - 1), vNsigma, kDetector::kTPC);
              break;
            }
          }
        } else {
          pidSelection = isPIDSelected(pidCut, std::vector<int>{3}, 4, confnSigmaAccetance.value, vNsigma, kDetector::kTPC);
        }
      } else {
        pidSelection = isPIDSelected(pidCut, std::vector<int>{3}, 4, confnSigmaAccetance.value, vNsigma, kDetector::kTPC);
      }
    } else {
      LOG(fatal) << "Other PID's are not supported by this trigger" << std::endl;
    }
    return pidSelection;
  }

  void init(o2::framework::InitContext&)
  {
    registry.add("fProcessedEvents", "CF Two Body - event filtered;;events}", HistType::kTH1F, {{2 + kLAST_CFTwoBodyTriggers, 0, 2 + kLAST_CFTwoBodyTriggers}});

    std::array<std::string, 2 + kLAST_CFTwoBodyTriggers> eventTitles = {"all", "rejected", "p-d", "l-d"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      registry.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    registry.add("fMultiplicityBefore", "Multiplicity before trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fMultiplicityAfter", "Multiplicity after trigger", HistType::kTH1F, {{1000, 0, 1000}});
    registry.add("fZvtxBefore", "Zvtx before trigger", HistType::kTH1F, {{1000, -15, 15}});
    registry.add("fZvtxAfter", "Zvtx after trigger", HistType::kTH1F, {{1000, -15, 15}});

    registry.add("fPtBeforeSel", "Transverse momentum of positive tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fEtaBeforeSel", "Pseudorapidity of positive tracks", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiBeforeSel", "Azimuthal angle of positive tracks", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

    registry.add("fPtAntiBeforeSel", "Transverse momentum of negative tracks", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fEtaAntiBeforeSel", "Pseudorapidity of negative tracks", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiAntiBeforeSel", "Azimuthal angle of negative tracks", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

    bool plotPerRadii = true;

    if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
      registry.add("fKstarPD", "CF - same event pd distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiPD", "CF - same event pd distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fEtaProtonAfterSel", "Pseudorapidity of Protons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
      registry.add("fPhiProtonAfterSel", "Azimuthal angle of Protons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

      registry.add("fPtAntiProtonAfterSel", "Transverse momentum of Protons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fEtaAntiProtonAfterSel", "Pseudorapidity of AntiProtons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
      registry.add("fPhiAntiProtonAfterSel", "Azimuthal angle of AntiProtons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

      closePairRejectionTT.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      registry.add("fKstarLD", "CF - same event ld distribution for particles;;events", HistType::kTH1F, {{8000, 0, 8}});
      registry.add("fKstarAntiLD", "CF - same event ld distribution for antiparticles;;events", HistType::kTH1F, {{8000, 0, 8}});

      registry.add("fPtLambdaAfterSel", "Transverse momentum of Lambdas which passed selections", HistType::kTH1F, {{1000, 0, 10}});
      registry.add("fPtAntiLambdaAfterSel", "Transverse momentum of AntidLambdas which passed selections", HistType::kTH1F, {{1000, 0, 10}});

      closePairRejectionTV0.init(&registry, &registryQA, ldeltaPhiMax, ldeltaEtaMax, plotPerRadii);
    }

    registry.add("fPtDeuteronAfterSel", "Transverse momentum of Deuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fEtaDeuteronAfterSel", "Pseudorapidity of Deuterons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiDeuteronAfterSel", "Azimuthal angle of Deuterons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});

    registry.add("fPtAntiDeuteronAfterSel", "Transverse momentum of Antideuterons which passed selections", HistType::kTH1F, {{1000, 0, 10}});
    registry.add("fEtaAntiDeuteronAfterSel", "Pseudorapidity of AntiDeuterons which passed selections", HistType::kTH1F, {{1000, -1, 1}});
    registry.add("fPhiAntiDeuteronAfterSel", "Azimuthal angle of AntiDeuterons which passed selections", HistType::kTH1F, {{1000, 0, TMath::TwoPi()}});
  }

  float mMassProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  float mMassDeuteron = o2::constants::physics::MassDeuteron;
  float mMassLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  void process(o2::aod::FemtoDreamCollision& col, o2::aod::FemtoDreamParticles& partsFemto)
  {
    // get partitions of all paritcles and antiparticles
    auto partsPD = partPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto partsAntiPD = partAntiPD->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    // get partions of V0s
    auto partsL = partL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());
    auto partsAntiL = partAntiL->sliceByCached(aod::femtodreamparticle::femtoDreamCollisionId, col.globalIndex());

    // magnetic field is need for close pair rejection
    auto magneticField = col.magField();

    registry.fill(HIST("fProcessedEvents"), 0);
    registry.fill(HIST("fMultiplicityBefore"), col.multV0M());
    registry.fill(HIST("fZvtxBefore"), col.posZ());

    // pass through the particles once to check if there are any particles of interest in the first place
    int Nproton = 0;
    int Nantiproton = 0;
    int Nlambda = 0;
    int Nantilambda = 0;
    int Ndeuteron = 0;
    int Nantideuteron = 0;

    for (auto pd : partsPD) {
      registry.fill(HIST("fPtBeforeSel"), pd.pt());
      registry.fill(HIST("fEtaBeforeSel"), pd.eta());
      registry.fill(HIST("fPhiBeforeSel"), pd.phi());

      // select protons
      if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
        if (SelectParticlePID(pd.pidcut(), o2::track::PID::Proton, pd.p(), pd.pt()) &&
            pd.pt() < confProtonPtMax.value &&
            pd.pt() > confProtonPtMin.value) {
          registry.fill(HIST("fPtProtonAfterSel"), pd.pt());
          registry.fill(HIST("fEtaProtonAfterSel"), pd.eta());
          registry.fill(HIST("fPhiProtonAfterSel"), pd.phi());
          Nproton++;
        }
      }
      // select deuterons
      if (SelectParticlePID(pd.pidcut(), o2::track::PID::Deuteron, pd.p(), pd.pt()) &&
          pd.pt() < confDeuteronPtMax.value &&
          pd.pt() > confDeuteronPtMin.value) {
        registry.fill(HIST("fPtDeuteronAfterSel"), pd.pt());
        registry.fill(HIST("fEtaDeuteronAfterSel"), pd.eta());
        registry.fill(HIST("fPhiDeuteronAfterSel"), pd.phi());
        Ndeuteron++;
      }
    }
    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      // select lambdas
      for (auto lambda : partsL) {
        registry.fill(HIST("fPtLambdaAfterSel"), lambda.pt());
        Nlambda++;
      }
    }

    for (auto antipd : partsAntiPD) {
      registry.fill(HIST("fPtAntiBeforeSel"), antipd.pt());
      registry.fill(HIST("fEtaAntiBeforeSel"), antipd.eta());
      registry.fill(HIST("fPhiAntiBeforeSel"), antipd.phi());
      // select antiprotons
      if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
        if (SelectParticlePID(antipd.pidcut(), o2::track::PID::Proton, antipd.p(), antipd.pt()) &&
            antipd.pt() < confProtonPtMax.value &&
            antipd.pt() > confProtonPtMin.value) {

          registry.fill(HIST("fPtAntiProtonAfterSel"), antipd.pt());
          registry.fill(HIST("fEtaAntiProtonAfterSel"), antipd.eta());
          registry.fill(HIST("fPhiAntiProtonAfterSel"), antipd.phi());
          Nantiproton++;
        }
      }
      // select antideuterons
      if (SelectParticlePID(antipd.pidcut(), o2::track::PID::Deuteron, antipd.p(), antipd.pt()) &&
          antipd.pt() < confDeuteronPtMax.value &&
          antipd.pt() > confDeuteronPtMin.value) {
        registry.fill(HIST("fPtAntiDeuteronAfterSel"), antipd.pt());
        registry.fill(HIST("fEtaAntiDeuteronAfterSel"), antipd.eta());
        registry.fill(HIST("fPhiAntiDeuteronAfterSel"), antipd.phi());
        Nantideuteron++;
      }
    }

    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      for (auto antilambda : partsAntiL) {
        // select antilambdas
        registry.fill(HIST("fPtAntiLambdaAfterSel"), antilambda.pt());
        Nlambda++;
      }
    }

    bool keepEvent[kLAST_CFTwoBodyTriggers] = {false, false};
    int lowKstarPairs[kLAST_CFTwoBodyTriggers] = {0, 0};

    bool pdPair = false;
    bool dpPair = false;
    double kStar = 0.;

    // trigger for pd pairs
    if (KstarTrigger.value == 0 || KstarTrigger.value == 11) {
      if (Ndeuteron > 0 && Nproton > 0) {
        // loop over all unique combinations of particles, excluding the self combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(partsPD, partsPD))) {

          // check if it is a pd pair
          // p1 => proton
          // p2 => deuteron
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Proton, p1.p(), p1.pt()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Deuteron, p2.p(), p1.pt()) &&
              p1.pt() < confProtonPtMax.value &&
              p1.pt() > confProtonPtMin.value &&
              p2.pt() < confDeuteronPtMax.value &&
              p2.pt() > confDeuteronPtMin.value) {
            pdPair = true;
          } else {
            pdPair = false;
          }

          // check if it is dp pair
          // p1 => deuteron
          // p2 => proton
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), p1.pt()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Proton, p2.p(), p2.pt()) &&
              p1.pt() < confDeuteronPtMax.value &&
              p1.pt() > confDeuteronPtMin.value &&
              p2.pt() < confProtonPtMax.value &&
              p2.pt() > confProtonPtMin.value) {
            dpPair = true;
          } else {
            dpPair = false;
          }

          // if neither is the case, skip
          if (!(pdPair || dpPair)) {
            continue;
          }

          // reject close pairs
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }

          // compute kstar depending on the pairing
          if (pdPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
          } else if (dpPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassDeuteron, p2, mMassProton);
          } else {
            kStar = confKstarTriggerLimit;
          }
          // check if the kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[kPD]++;
            registry.fill(HIST("fKstarPD"), kStar);
          }
        }
      }

      if (Nantideuteron > 0 && Nantiproton > 0) {
        // loop over all unique combinations of antiparticles, excluding the self combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(partsAntiPD, partsAntiPD))) {

          // check if it is a (anti)pd pair
          // p1 => antiproton
          // p2 => antideuteron
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Proton, p1.p(), p1.pt()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Deuteron, p2.p(), p2.pt()) &&
              p1.pt() < confProtonPtMax.value &&
              p1.pt() > confProtonPtMin.value &&
              p2.pt() < confDeuteronPtMax.value &&
              p2.pt() > confDeuteronPtMin.value) {
            pdPair = true;
          } else {
            pdPair = false;
          }

          // check if it is (anti)dp pair
          // p1 => antideuteron
          // p2 => antiproton
          if (SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), p1.pt()) &&
              SelectParticlePID(p2.pidcut(), o2::track::PID::Proton, p2.p(), p2.pt()) &&
              p1.pt() < confDeuteronPtMax.value &&
              p1.pt() > confDeuteronPtMin.value &&
              p2.pt() < confProtonPtMax.value &&
              p2.pt() > confProtonPtMin.value) {
            dpPair = true;
          } else {
            dpPair = false;
          }

          // if neither is the case, skip
          if (!(pdPair || dpPair)) {
            continue;
          }

          // reject close pairs
          if (closePairRejectionTT.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }

          // compute kstar depending on the pairing
          if (pdPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassProton, p2, mMassDeuteron);
          } else if (dpPair) {
            kStar = FemtoDreamMath::getkstar(p1, mMassDeuteron, p2, mMassProton);
          } else {
            kStar = confKstarTriggerLimit;
          }

          // check if the kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[kPD]++;
            registry.fill(HIST("fKstarAntiPD"), kStar);
          }
        }
      }
    }

    // trigger for ld pairs
    if (KstarTrigger.value == 1 || KstarTrigger.value == 11) {
      if (Ndeuteron > 0 && Nlambda > 0) {
        // loop over all unique combinations
        for (auto& [p1, p2] : combinations(soa::CombinationsUpperIndexPolicy(partsPD, partsL))) {
          // check if the particle is a deuteron
          // we do not need to check the V0s
          if (!SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), p1.pt())) {
            continue;
          }
          if (closePairRejectionTV0.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          kStar = FemtoDreamMath::getkstar(p1, mMassLambda, p2, mMassDeuteron);
          // check if kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[1]++;
            registry.fill(HIST("fKstarLD"), kStar);
          }
        }
      }
      if (Nantideuteron > 0 && Nantilambda > 0) {
        for (auto& [p1, p2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(partsAntiPD, partsAntiL))) {
          // check if the particle is a antideuteron
          if (!SelectParticlePID(p1.pidcut(), o2::track::PID::Deuteron, p1.p(), p1.pt())) {
            continue;
          }
          if (closePairRejectionTV0.isClosePair(p1, p2, partsFemto, magneticField)) {
            continue;
          }
          auto kstar = FemtoDreamMath::getkstar(p1, mMassLambda, p2, mMassDeuteron);
          // check if kstar is below threshold
          if (kStar < confKstarTriggerLimit.value) {
            lowKstarPairs[1]++;
            registry.fill(HIST("fKstarAntiLD"), kstar);
          }
        }
      }
    }

    // if we found any pair below the kstar limit, keep the event
    if (lowKstarPairs[kPD] > 0) {
      keepEvent[kPD] = true;
      registry.fill(HIST("fProcessedEvents"), 2);
    }
    if (lowKstarPairs[kLD] > 0) {
      keepEvent[kLD] = true;
      registry.fill(HIST("fProcessedEvents"), 3);
    }

    // fill table for the trigger
    tags(keepEvent[kPD], keepEvent[kLD]);

    if (keepEvent[kPD] > 0 || keepEvent[kLD] > 0) {
      registry.fill(HIST("fMultiplicityAfter"), col.multV0M());
      registry.fill(HIST("fZvtxAfter"), col.posZ());
    } else {
      registry.fill(HIST("fProcessedEvents"), 1);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<CFFilterTwoN>(cfg)};
}
