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

/// \file TrackQa.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for tracks
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
#include "PWGCF/FemtoUnited/Core/VzeroHistManager.h"

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

  struct : ConfigurableGroup {
    std::string prefix = std::string("CollisionSelection");
    Configurable<float> vtxZMin{"vtxZMin", -10., "Minimum vertex Z position (cm)"};
    Configurable<float> vtxZMax{"vtxZMax", 10., "Maximum vertex Z position (cm)"};
    Configurable<float> multMin{"multMin", 0, "Minimum multiplicity"};
    Configurable<float> multMax{"multMax", 200, "Maximum multiplicity"};
    Configurable<float> centMin{"centMin", 0.0f, "Minimum centrality (multiplicity percentile)"};
    Configurable<float> centMax{"centMax", 100.0f, "Maximum centrality (multiplicity percentile)"};
    Configurable<float> spherMin{"spherMin", 0.0f, "Minimum centrality (multiplicity percentile)"};
    Configurable<float> spherMax{"spherMax", 2.0f, "Maximum centrality (multiplicity percentile)"};
    Configurable<float> magFieldMin{"magFieldMin", -1.0f, "Minimum magnetic field strength (T)"};
    Configurable<float> magFieldMax{"magFieldMax", 1.0f, "Maximum magnetic field strength (T)"};
  } CollisionSelection;

  Filter filterVtxz = femtocollisions::posZ >= CollisionSelection.vtxZMin && femtocollisions::posZ <= CollisionSelection.vtxZMax;
  Filter filterMult = femtocollisions::mult >= CollisionSelection.multMin && femtocollisions::mult <= CollisionSelection.multMax;
  Filter filterCent = femtocollisions::cent >= CollisionSelection.centMin && femtocollisions::cent <= CollisionSelection.centMax;
  Filter filterSpher = femtocollisions::sphericity >= CollisionSelection.centMin && femtocollisions::sphericity <= CollisionSelection.centMax;
  Filter filterMagField = femtocollisions::magField >= CollisionSelection.magFieldMin && femtocollisions::magField <= CollisionSelection.magFieldMax;

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  struct : ConfigurableGroup {
    std::string prefix = std::string("CollisionBinning");
    ConfigurableAxis vtZ{"vtZ", {200, -10, 10}, "Vertex Z binning"};
    ConfigurableAxis mult{"mult", {200, 0, 200}, "Multiplicity binning"};
    ConfigurableAxis cent{"cent", {100, 0.0f, 100.0f}, "Centrality (multiplicity percentile) binning"};
    ConfigurableAxis spher{"spher", {200, 0.0f, 2.0f}, "Sphericity binning"};
    ConfigurableAxis magField{"magField", {2, -1, 1}, "Magnetic field binning"};
  } CollisionBinning;

  struct : ConfigurableGroup {
    std::string prefix = std::string("VzeroSelection");
    Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
    Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
    Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
    Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
    Configurable<float> phiMin{"phiMin", 0.f, "Minimum eta"};
    Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
    Configurable<float> massMin{"massMin", 1.f, "Minimum invariant mass"};
    Configurable<float> massMax{"massMax", 1.2f, "Maximum invariant mass"};
    Configurable<femtodatatypes::VzeroMaskType> mask{"mask", 1, "Bitmask for V0 selection"};
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

  SliceCache cache;

  using V0s = o2::soa::Join<FUVzeros, FUVzeroMasks, FUVzeroDaus, FUVzeroExtras, FUVzeroDauExts>;

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
    ncheckbit(femtovzeros::vzeroMask, VzeroSelection.mask);

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
      {colhistmanager::kPosz, {CollisionBinning.vtZ}},
      {colhistmanager::kMult, {CollisionBinning.mult}},
      {colhistmanager::kCent, {CollisionBinning.cent}},
      {colhistmanager::kSphericity, {CollisionBinning.spher}},
      {colhistmanager::kMagField, {CollisionBinning.magField}},
      {colhistmanager::kPoszVsMult, {CollisionBinning.vtZ, CollisionBinning.mult}},
      {colhistmanager::kPoszVsCent, {CollisionBinning.vtZ, CollisionBinning.cent}},
      {colhistmanager::kCentVsMult, {CollisionBinning.cent, CollisionBinning.mult}},
      {colhistmanager::kMultVsSphericity, {CollisionBinning.mult, CollisionBinning.spher}},
      {colhistmanager::kCentVsSphericity, {CollisionBinning.cent, CollisionBinning.spher}}};
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
      {vzerohistmanager::kSign, {VzeroBinning.sign}},
      {vzerohistmanager::kDecayDauDca, {VzeroBinning.dauDca}},
      {vzerohistmanager::kDecayVtxX, {VzeroBinning.decayVertex}},
      {vzerohistmanager::kDecayVtxY, {VzeroBinning.decayVertex}},
      {vzerohistmanager::kDecayVtxZ, {VzeroBinning.decayVertex}},
      {vzerohistmanager::kTransRadius, {VzeroBinning.transRadius}},
      {vzerohistmanager::kKaonMass, {VzeroBinning.kaonMass}},
      {vzerohistmanager::kPtVsEta, {VzeroBinning.pt, VzeroBinning.eta}},
      {vzerohistmanager::kPtVsPhi, {VzeroBinning.pt, VzeroBinning.phi}},
      {vzerohistmanager::kPhiVsEta, {VzeroBinning.phi, VzeroBinning.eta}},
      {vzerohistmanager::kPosDaughTpcCluster, {VzeroBinning.dauTpcCluster}},
      {vzerohistmanager::kPosDaughPtVsDcaxy, {VzeroBinning.dauPt, VzeroBinning.dauDcaxy}},
      {vzerohistmanager::kPosDaughPtVsDcaz, {VzeroBinning.dauPt, VzeroBinning.dauDcaz}},
      {vzerohistmanager::kPosDaughPtVsDca, {VzeroBinning.dauPt, VzeroBinning.dauDca}},
      {vzerohistmanager::kPosDaughTpcNsigma, {VzeroBinning.dauP, VzeroBinning.dauTpcNsigma}},
      {vzerohistmanager::kNegDaughTpcCluster, {VzeroBinning.dauTpcCluster}},
      {vzerohistmanager::kNegDaughPtVsDcaxy, {VzeroBinning.dauPt, VzeroBinning.dauDcaxy}},
      {vzerohistmanager::kNegDaughPtVsDcaz, {VzeroBinning.dauPt, VzeroBinning.dauDcaz}},
      {vzerohistmanager::kNegDaughPtVsDca, {VzeroBinning.dauPt, VzeroBinning.dauDca}},
      {vzerohistmanager::kNegDaughTpcNsigma, {VzeroBinning.dauP, VzeroBinning.dauTpcNsigma}},
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
