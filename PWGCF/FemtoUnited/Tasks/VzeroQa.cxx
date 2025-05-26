// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file VzeroQa.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for vzeros
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Expressions.h"
#include "CommonConstants/MathConstants.h"

#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoVzerosDerived.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/VzeroHistManager.h"
#include "PWGCF/FemtoUnited/Core/VzeroSelection.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

struct VzeroQa {

  struct : ConfigurableGroup {
    std::string prefix = std::string("Options");
    Configurable<bool> correlatedPlots{"correlatedPlots", false, "Enable multidimensional histogramms. High memory consumption."};
  } Options;

  collisionselection::ConfCollisionSelection collisionSelection;
  Filter filterVtxz = femtocollisions::posZ >= collisionSelection.vtxZMin && femtocollisions::posZ <= collisionSelection.vtxZMax;
  Filter filterMult = femtocollisions::mult >= collisionSelection.multMin && femtocollisions::mult <= collisionSelection.multMax;
  Filter filterCent = femtocollisions::cent >= collisionSelection.centMin && femtocollisions::cent <= collisionSelection.centMax;
  Filter filterSpher = femtocollisions::sphericity >= collisionSelection.centMin && femtocollisions::sphericity <= collisionSelection.centMax;
  Filter filterMagField = femtocollisions::magField >= collisionSelection.magFieldMin && femtocollisions::magField <= collisionSelection.magFieldMax;

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  colhistmanager::ConfCollisionBinning collisionBinning;

  using V0s = o2::soa::Join<FUVzeros, FUVzeroMasks, FUVzeroDaus, FUVzeroExtras, FUVzeroDauExts>;

  SliceCache cache;

  struct : ConfigurableGroup {
    std::string prefix = std::string("VzeroSelection");
    Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
    Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
    Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
    Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
    Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
    Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
    Configurable<float> massMin{"massMin", 1.f, "Minimum invariant mass for Lambda"};
    Configurable<float> massMax{"massMax", 1.2f, "Maximum invariant mass for Lambda"};
    Configurable<float> antiMassMin{"antiMassMin", 0, "Minimum invariant mass for AntiLambda"};
    Configurable<float> antiMassMax{"antiMassMax", 9999, "Maximum invariant mass for AntiLambda"};
    Configurable<o2::aod::femtodatatypes::VzeroMaskType> mask{"mask", 0, "Bitmask for V0 selection"};
    Configurable<o2::aod::femtodatatypes::VzeroDauPidMaskType> dauPidMask{"dauPidMask", 6, "Bitmask for V0 selection"};
  } VzeroSelection;

  struct : ConfigurableGroup {
    std::string prefix = std::string("VzeroDaughterSelection");
    Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
    Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
    Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
    Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
    Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
    Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  } VzeroDaughterSelection;

  Partition<V0s> VzeroPartition =
    (femtobase::pt > VzeroSelection.ptMin) &&
    (femtobase::pt < VzeroSelection.ptMax) &&
    (femtobase::eta > VzeroSelection.etaMin) &&
    (femtobase::eta < VzeroSelection.etaMax) &&
    (femtobase::phi > VzeroSelection.phiMin) &&
    (femtobase::phi < VzeroSelection.phiMax) &&
    (femtovzeros::vzeroMass > VzeroSelection.massMin) &&
    (femtovzeros::vzeroMass < VzeroSelection.massMax) &&
    (femtovzeros::posDauPt > VzeroDaughterSelection.ptMin) &&
    (femtovzeros::posDauPt < VzeroDaughterSelection.ptMax) &&
    (femtovzeros::posDauEta > VzeroDaughterSelection.etaMin) &&
    (femtovzeros::posDauEta < VzeroDaughterSelection.etaMax) &&
    (femtovzeros::posDauPhi > VzeroDaughterSelection.phiMin) &&
    (femtovzeros::posDauPhi < VzeroDaughterSelection.phiMax) &&
    (femtovzeros::negDauPt > VzeroDaughterSelection.ptMin) &&
    (femtovzeros::negDauPt < VzeroDaughterSelection.ptMax) &&
    (femtovzeros::negDauEta > VzeroDaughterSelection.etaMin) &&
    (femtovzeros::negDauEta < VzeroDaughterSelection.etaMax) &&
    (femtovzeros::negDauPhi > VzeroDaughterSelection.phiMin) &&
    (femtovzeros::negDauPhi < VzeroDaughterSelection.phiMax) &&
    ncheckbit(femtovzeros::vzeroMask, VzeroSelection.mask) &&
    ncheckbit(femtovzeros::vzeroDauPidMask, VzeroSelection.dauPidMask);
  Preslice<V0s> perColReco = aod::femtobase::collisionId;

  struct : ConfigurableGroup {
    std::string prefix = std::string("VzeroBinning");
    ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
    ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
    ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
    ConfigurableAxis mass{"mass", {{200, 1.0, 1.2}}, "Mass"};
    ConfigurableAxis dauPt{"dauPt", {{600, 0, 6}}, "daughter pt"};
    ConfigurableAxis dauEta{"dauEta", {{300, -1.5, 1.5}}, "daugher eta"};
    ConfigurableAxis dauPhi{"dauPhi", {{720, 0., 1.f * o2::constants::math::TwoPI}}, "daughter phi"};
    ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
    ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
    ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
    ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
    ConfigurableAxis kaonMass{"kaonMass", {{100, 0.45, 0.55}}, "Mass for kaon hypothesis"};
    ConfigurableAxis dauTpcCluster{"dauTpcCluster", {{153, -0.5, 152.5}}, "TPC cluster of daughters"};
    ConfigurableAxis dauP{"dauP", {{600, 0, 6}}, "Momentum binning for TPC Nsigma of daughters"};
    ConfigurableAxis dauDcaxy{"dauDcaxy", {{300, -1.5, 1.5}}, "DCAxy for daughters"};
    ConfigurableAxis dauDcaz{"dauDcaz", {{300, -1.5, 1.5}}, "Dcaz for daughters"};
    ConfigurableAxis dauDca{"dauDca", {{150, 0, 0.3}}, "Dca for daughters"};
    ConfigurableAxis dauTpcNsigma{"dauTpcNsigma", {{600, -6, 6}}, "TPC Nsigma for daughters"};
  } VzeroBinning;

  HistogramRegistry hRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  colhistmanager::CollisionHistManager colHistManager;
  vzerohistmanager::VzeroHistManager vzeroHistManager;

  void init(InitContext&)
  {
    // create a map for histogram specs
    std::map<colhistmanager::ColHist, std::vector<framework::AxisSpec>> colHistSpec = {
      {colhistmanager::kPosz, {collisionBinning.vtZ}},
      {colhistmanager::kMult, {collisionBinning.mult}},
      {colhistmanager::kCent, {collisionBinning.cent}},
      {colhistmanager::kSphericity, {collisionBinning.spher}},
      {colhistmanager::kMagField, {collisionBinning.magField}},
      {colhistmanager::kPoszVsMult, {collisionBinning.vtZ, collisionBinning.mult}},
      {colhistmanager::kPoszVsCent, {collisionBinning.vtZ, collisionBinning.cent}},
      {colhistmanager::kCentVsMult, {collisionBinning.cent, collisionBinning.mult}},
      {colhistmanager::kMultVsSphericity, {collisionBinning.mult, collisionBinning.spher}},
      {colhistmanager::kCentVsSphericity, {collisionBinning.cent, collisionBinning.spher}}};
    colHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, colHistSpec);

    std::map<vzerohistmanager::VzeroHist, std::vector<framework::AxisSpec>> vzeroHistSpec = {
      {vzerohistmanager::kPt, {VzeroBinning.pt}},
      {vzerohistmanager::kEta, {VzeroBinning.eta}},
      {vzerohistmanager::kPhi, {VzeroBinning.phi}},
      {vzerohistmanager::kMass, {VzeroBinning.mass}},
      {vzerohistmanager::kPosDauPt, {VzeroBinning.dauPt}},
      {vzerohistmanager::kPosDauEta, {VzeroBinning.dauEta}},
      {vzerohistmanager::kPosDauPhi, {VzeroBinning.dauPhi}},
      {vzerohistmanager::kNegDauPt, {VzeroBinning.dauPt}},
      {vzerohistmanager::kNegDauEta, {VzeroBinning.dauEta}},
      {vzerohistmanager::kNegDauPhi, {VzeroBinning.dauPhi}},
      {vzerohistmanager::kDecayDauDca, {VzeroBinning.dauDca}},
      {vzerohistmanager::kDecayVtxX, {VzeroBinning.decayVertex}},
      {vzerohistmanager::kDecayVtxY, {VzeroBinning.decayVertex}},
      {vzerohistmanager::kDecayVtxZ, {VzeroBinning.decayVertex}},
      {vzerohistmanager::kTransRadius, {VzeroBinning.transRadius}},
      {vzerohistmanager::kKaonMass, {VzeroBinning.kaonMass}},
      {vzerohistmanager::kPtVsEta, {VzeroBinning.pt, VzeroBinning.eta}},
      {vzerohistmanager::kPtVsPhi, {VzeroBinning.pt, VzeroBinning.phi}},
      {vzerohistmanager::kPhiVsEta, {VzeroBinning.phi, VzeroBinning.eta}},
      {vzerohistmanager::kPosDauTpcCluster, {VzeroBinning.dauTpcCluster}},
      {vzerohistmanager::kPosDauPtVsDcaxy, {VzeroBinning.dauPt, VzeroBinning.dauDcaxy}},
      {vzerohistmanager::kPosDauPtVsDcaz, {VzeroBinning.dauPt, VzeroBinning.dauDcaz}},
      {vzerohistmanager::kPosDauPtVsDca, {VzeroBinning.dauPt, VzeroBinning.dauDca}},
      {vzerohistmanager::kPosDauProtonTpcNsigma, {VzeroBinning.dauP, VzeroBinning.dauTpcNsigma}},
      {vzerohistmanager::kPosDauPionTpcNsigma, {VzeroBinning.dauP, VzeroBinning.dauTpcNsigma}},
      {vzerohistmanager::kNegDauTpcCluster, {VzeroBinning.dauTpcCluster}},
      {vzerohistmanager::kNegDauPtVsDcaxy, {VzeroBinning.dauPt, VzeroBinning.dauDcaxy}},
      {vzerohistmanager::kNegDauPtVsDcaz, {VzeroBinning.dauPt, VzeroBinning.dauDcaz}},
      {vzerohistmanager::kNegDauPtVsDca, {VzeroBinning.dauPt, VzeroBinning.dauDca}},
      {vzerohistmanager::kNegDauProtonTpcNsigma, {VzeroBinning.dauP, VzeroBinning.dauTpcNsigma}},
      {vzerohistmanager::kNegDauPionTpcNsigma, {VzeroBinning.dauP, VzeroBinning.dauTpcNsigma}},
    };

    vzeroHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, vzeroHistSpec);
  };

  void process(FilteredCollision const& col, V0s const& /*v0s*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto vzeroSlice = VzeroPartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    for (auto const& vzero : vzeroSlice) {
      vzeroHistManager.fill<modes::Mode::kANALYSIS_QA>(vzero);
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<VzeroQa>(cfgc),
  };
  return workflow;
}
