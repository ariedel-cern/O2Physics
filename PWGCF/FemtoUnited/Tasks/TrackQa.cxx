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
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/TrackHistManager.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

struct TrackQa {

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
    std::string prefix = std::string("TrackSelection");
    Configurable<int> pdgCode{"pdgCode", 2212, "Track PDG code"};
    Configurable<float> ptMin{"ptMin", 0.f, "Minimum pT"};
    Configurable<float> ptMax{"ptMax", 999.f, "Maximum pT"};
    Configurable<float> etaMin{"etaMin", -10.f, "Minimum eta"};
    Configurable<float> etaMax{"etaMax", 10.f, "Maximum eta"};
    Configurable<femtodatatypes::TrackMaskType> mask{"mask", 0, "Bitmask for track selection"};
    Configurable<femtodatatypes::TrackTPCMaskType> tpcMask{"tpcMask", 0, "Bitmask for TPC PID selection"};
    Configurable<femtodatatypes::TrackTOFMaskType> tofMask{"tofMask", 0, "Bitmask for TOF PID selection"};
    Configurable<femtodatatypes::TrackTPCTOFMaskType> tpctofMask{"tpctofMask", 0, "Bitmask for TPC+TOF selection"};
    Configurable<float> pidThres{"pidThres", 10.f, "Momentum threshold for PID of tracks with large momentum"};
    Configurable<bool> pidSwitch{"pidSwitch", true, "IF switch is true, use TPC+TOF PID for large momentum tracks and use only TOF PID when set to false"};
  } TrackSelections;

  SliceCache cache;

  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  Partition<Tracks> TrackPartition =
    (femtobase::pt > TrackSelections.ptMin) &&
    (femtobase::pt < TrackSelections.ptMax) &&
    (femtobase::eta > TrackSelections.etaMin) &&
    (femtobase::eta < TrackSelections.etaMax) &&
    ncheckbit(femtotracks::trackMask, TrackSelections.mask) &&
    ifnode(femtobase::pt * (nexp(femtobase::eta) + nexp(-1.f * femtobase::eta)) / 2.f <= TrackSelections.pidThres,
           ncheckbit(femtotracks::tpcMask, TrackSelections.tpcMask),
           ifnode(TrackSelections.pidSwitch, ncheckbit(femtotracks::tpctofMask, TrackSelections.tpctofMask), ncheckbit(femtotracks::tofMask, TrackSelections.tofMask)));

  Preslice<Tracks> perColReco = aod::femtobase::collisionId;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackBinning");
    ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
    ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
    ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
    ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
    ConfigurableAxis itsCluster{"itsCluster", {{8, -0.5, 7.5}}, "ITS cluster"};
    ConfigurableAxis itsClusterIb{"itsClusterIb", {{4, -0.5, 3.5}}, "ITS cluster in inner barrel"};
    ConfigurableAxis tpcCluster{"tpcCluster", {{153, -0.5, 152.5}}, "TPC cluster"};
    ConfigurableAxis tpcClusterShared{"tpcClusterShared", {{153, -0.5, 152.5}}, "TPC cluster shared"};
  } TrackBinning;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackPidBinning");
    ConfigurableAxis p{"p", {{300, 0, 6}}, "Momentum axis"};
    ConfigurableAxis tpcSignal{"tpcSignal", {{150, 0, 150}}, "TPC Signal"};
    ConfigurableAxis tpcElectron{"tpcElectron", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tpcPion{"tpcPion", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tpcKaon{"tpcKaon", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tpcProton{"tpcProton", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tpcDeuteron{"tpcDeuteron", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tpcTriton{"tpcTriton", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tpcHelium{"tpcHelium", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tofBeta{"tofBeta", {{150, 0, 1.5}}, "TPC Signal"};
    ConfigurableAxis tofElectron{"tofElectron", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tofPion{"tofPion", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tofKaon{"tofKaon", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tofProton{"tofProton", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tofDeuteron{"tofDeuteron", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tofTriton{"tofTriton", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tofHelium{"tofHelium", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tpctofElectron{"tpctofElectron", {{300, 0, 3}}, "tpctof PID for proton"};
    ConfigurableAxis tpctofPion{"tpctofPion", {{300, 0, 3}}, "tpctof PID for proton"};
    ConfigurableAxis tpctofKaon{"tpctofKaon", {{300, 0, 3}}, "tpctof PID for proton"};
    ConfigurableAxis tpctofProton{"tpctofProton", {{300, 0, 3}}, "tpctof PID for proton"};
    ConfigurableAxis tpctofDeuteron{"tpctofDeuteron", {{300, 0, 3}}, "tpctof PID for proton"};
    ConfigurableAxis tpctofTriton{"tpctofTriton", {{300, 0, 3}}, "tpctof PID for proton"};
    ConfigurableAxis tpctofHelium{"tpctofHelium", {{300, 0, 3}}, "tpctof PID for proton"};
  } TrackPidBinning;

  HistogramRegistry hRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  colhistmanager::CollisionHistManager colHistManager;
  trackhistmanager::TrackHistManager trackHistManager;

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

    std::map<trackhistmanager::TrackHist, std::vector<framework::AxisSpec>> trackHistSpec = {
      {trackhistmanager::kPt, {TrackBinning.pt}},
      {trackhistmanager::kEta, {TrackBinning.eta}},
      {trackhistmanager::kPhi, {TrackBinning.phi}},
      {trackhistmanager::kSign, {TrackBinning.sign}},
      {trackhistmanager::kItsCluster, {TrackBinning.itsCluster}},
      {trackhistmanager::kItsClusterIb, {TrackBinning.itsClusterIb}},
      {trackhistmanager::kPtVsEta, {TrackBinning.pt, TrackBinning.eta}},
      {trackhistmanager::kPtVsPhi, {TrackBinning.pt, TrackBinning.phi}},
      {trackhistmanager::kPhiVsEta, {TrackBinning.phi, TrackBinning.eta}},
      {trackhistmanager::kPtVsItsCluster, {TrackBinning.pt, TrackBinning.itsCluster}},
      {trackhistmanager::kPtVsTpcCluster, {TrackBinning.pt, TrackBinning.tpcCluster}},
      {trackhistmanager::kPtVsTpcClusterShared, {TrackBinning.pt, TrackBinning.tpcClusterShared}},
      {trackhistmanager::kTpcClusterVsTpcClusterShared, {TrackBinning.tpcCluster, TrackBinning.tpcClusterShared}},
      {trackhistmanager::kTpcCluster, {TrackBinning.tpcCluster}},
      {trackhistmanager::kTpcClusterShared, {TrackBinning.tpcClusterShared}},
      {trackhistmanager::kTpcSignal, {TrackPidBinning.p, TrackPidBinning.tpcSignal}},
      {trackhistmanager::kTpcElectron, {TrackPidBinning.p, TrackPidBinning.tpcElectron}},
      {trackhistmanager::kTpcPion, {TrackPidBinning.p, TrackPidBinning.tpcPion}},
      {trackhistmanager::kTpcKaon, {TrackPidBinning.p, TrackPidBinning.tpcKaon}},
      {trackhistmanager::kTpcProton, {TrackPidBinning.p, TrackPidBinning.tpcProton}},
      {trackhistmanager::kTpcDeuteron, {TrackPidBinning.p, TrackPidBinning.tpcDeuteron}},
      {trackhistmanager::kTpcTriton, {TrackPidBinning.p, TrackPidBinning.tpcTriton}},
      {trackhistmanager::kTpcHelium, {TrackPidBinning.p, TrackPidBinning.tpcHelium}},
      {trackhistmanager::kTofBeta, {TrackPidBinning.p, TrackPidBinning.tofBeta}},
      {trackhistmanager::kTofElectron, {TrackPidBinning.p, TrackPidBinning.tofElectron}},
      {trackhistmanager::kTofPion, {TrackPidBinning.p, TrackPidBinning.tofPion}},
      {trackhistmanager::kTofKaon, {TrackPidBinning.p, TrackPidBinning.tofKaon}},
      {trackhistmanager::kTofProton, {TrackPidBinning.p, TrackPidBinning.tofProton}},
      {trackhistmanager::kTofDeuteron, {TrackPidBinning.p, TrackPidBinning.tofDeuteron}},
      {trackhistmanager::kTofTriton, {TrackPidBinning.p, TrackPidBinning.tofTriton}},
      {trackhistmanager::kTofHelium, {TrackPidBinning.p, TrackPidBinning.tofHelium}},
      {trackhistmanager::kTpctofElectron, {TrackPidBinning.p, TrackPidBinning.tpctofElectron}},
      {trackhistmanager::kTpctofPion, {TrackPidBinning.p, TrackPidBinning.tpctofPion}},
      {trackhistmanager::kTpctofKaon, {TrackPidBinning.p, TrackPidBinning.tpctofKaon}},
      {trackhistmanager::kTpctofProton, {TrackPidBinning.p, TrackPidBinning.tpctofProton}},
      {trackhistmanager::kTpctofDeuteron, {TrackPidBinning.p, TrackPidBinning.tpctofDeuteron}},
      {trackhistmanager::kTpctofTriton, {TrackPidBinning.p, TrackPidBinning.tpctofTriton}},
      {trackhistmanager::kTpctofHelium, {TrackPidBinning.p, TrackPidBinning.tpctofHelium}}};

    trackHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, trackHistSpec);
  };

  void process(FilteredCollision const& col, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto trackSlice = TrackPartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    for (auto const& track : trackSlice) {
      trackHistManager.fill<modes::Mode::kANALYSIS_QA>(track);
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TrackQa>(cfgc),
  };
  return workflow;
}
