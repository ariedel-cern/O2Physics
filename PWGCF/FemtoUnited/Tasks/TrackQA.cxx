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

/// \file TrackQA.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Expressions.h"

#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/TrackHistManager.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct femtounitedtrackqa {

  struct : ConfigurableGroup {
    std::string prefix = std::string("Options");
    Configurable<bool> CorrelatedPlots{"CorrelatedPlots", false, "Enable multidimensional histogramms. High memory consumption."};
  } Options;

  struct : ConfigurableGroup {
    std::string prefix = std::string("CollisionSelection");
    Configurable<float> VtxZMin{"VtxZMin", -10., "Minimum vertex Z position (cm)"};
    Configurable<float> VtxZMax{"VtxZMax", 10., "Maximum vertex Z position (cm)"};
    Configurable<float> MultMin{"MultMin", 0, "Minimum multiplicity"};
    Configurable<float> MultMax{"MultMax", 200, "Maximum multiplicity"};
    Configurable<float> CentMin{"CentMin", 0.0f, "Minimum centrality (multiplicity percentile)"};
    Configurable<float> CentMax{"CentMax", 100.0f, "Maximum centrality (multiplicity percentile)"};
    Configurable<float> MagFieldMin{"MagFieldMin", -1.0f, "Minimum magnetic field strength (T)"};
    Configurable<float> MagFieldMax{"MagFieldMax", 1.0f, "Maximum magnetic field strength (T)"};
  } CollisionSelection;

  Filter VtxZ = femtocollisions::posZ >= CollisionSelection.VtxZMin && femtocollisions::posZ <= CollisionSelection.VtxZMax;
  Filter Mult = femtocollisions::mult >= CollisionSelection.MultMin && femtocollisions::mult <= CollisionSelection.MultMax;
  Filter Cent = femtocollisions::cent >= CollisionSelection.CentMin && femtocollisions::cent <= CollisionSelection.CentMax;
  Filter MagField = femtocollisions::magField >= CollisionSelection.MagFieldMin && femtocollisions::magField <= CollisionSelection.MagFieldMax;

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  struct : ConfigurableGroup {
    std::string prefix = std::string("CollisionBinning");
    ConfigurableAxis VtZ{"VtZ", {200, -10, 10}, "Vertex Z binning"};
    ConfigurableAxis Mult{"Mult", {200, 0, 200}, "Multiplicity binning"};
    ConfigurableAxis Cent{"Cent", {100, 0.0f, 100.0f}, "Centrality (multiplicity percentile) binning"};
    ConfigurableAxis MagField{"MagField", {2, -1, 1}, "Magnetic field binning"};
  } CollisionBinning;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackSelection");
    Configurable<int> PDGCode{"PDGCode", 2212, "Track PDG code"};
    Configurable<float> PtMin{"PtMin", 0.f, "Minimum pT"};
    Configurable<float> PtMax{"PtMax", 999.f, "Maximum pT"};
    Configurable<float> EtaMin{"EtaMin", -10.f, "Minimum eta"};
    Configurable<float> EtaMax{"EtaMax", 10.f, "Maximum eta"};
    Configurable<femtodatatypes::TrackMaskType> Mask{"TrackMask", 0, "Bitmask for track selection"};
    Configurable<femtodatatypes::TrackTPCMaskType> TPCMask{"TrackTPCMask", 0, "Bitmask for TPC PID selection"};
    Configurable<femtodatatypes::TrackTOFMaskType> TOFMask{"TrackTOFMask", 0, "Bitmask for TOF PID selection"};
    Configurable<femtodatatypes::TrackTPCTOFMaskType> TPCTOFMask{"TrackTPCTOFMask", 0, "Bitmask for TPC+TOF selection"};
    Configurable<float> PIDThres{"PIDThres", 10.f, "PID Threshold"};
  } TrackSelections;

  SliceCache cache;

  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks>; // , FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  Partition<Tracks> TrackPartition =
    (femtobase::pt > TrackSelections.PtMin) &&
    (femtobase::pt < TrackSelections.PtMax) &&
    (femtobase::eta > TrackSelections.EtaMin) &&
    (femtobase::eta < TrackSelections.EtaMax) &&
    ncheckbit(femtotracks::trackMask, TrackSelections.Mask) &&
    ifnode(femtobase::pt * (nexp(femtobase::eta) + nexp(-1.f * femtobase::eta)) / 2.f <= TrackSelections.PIDThres,
           ncheckbit(femtotracks::tpcMask, TrackSelections.TPCMask),
           ncheckbit(femtotracks::tpctofMask, TrackSelections.TPCTOFMask));

  Preslice<Tracks> perColReco = aod::femtobase::collisionId;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackBinning");
    ConfigurableAxis Pt{"Pt", {{600, 0, 6}}, "Pt binning"};
    ConfigurableAxis Eta{"Eta", {{300, -1.5, 1.5}}, "Eta binning"};
    ConfigurableAxis Phi{"Phi", {{720, 0, TMath::TwoPi()}}, "Phi binning"};
  } TrackBinning;

  HistogramRegistry HRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  femtounitedcolhistmanager::CollisionHistManager ColHists;
  femtounitedtrackhistmanager::TrackHistManager TrackHists;

  void init(InitContext&)
  {
    // create a map for histogram specs
    std::map<femtounitedcolhistmanager::EventVariable, std::vector<framework::AxisSpec>> CollisionHistogramSpec = {
      {femtounitedcolhistmanager::kPosZ, {CollisionBinning.VtZ}},
      {femtounitedcolhistmanager::kMult, {CollisionBinning.Mult}},
      {femtounitedcolhistmanager::kCent, {CollisionBinning.Cent}},
      {femtounitedcolhistmanager::kMagField, {CollisionBinning.MagField}}};
    ColHists.init<femtounitedcolhistmanager::Mode::kANALYSIS>(&HRegistry, CollisionHistogramSpec);

    std::map<femtounitedtrackhistmanager::TrackVariable, std::vector<framework::AxisSpec>> TrackHistogramSpec = {
      {femtounitedtrackhistmanager::kPt, {TrackBinning.Pt}},
      {femtounitedtrackhistmanager::kEta, {TrackBinning.Eta}},
      {femtounitedtrackhistmanager::kPhi, {TrackBinning.Phi}}};
    TrackHists.init<femtounitedtrackhistmanager::Mode::kANALYSIS>(&HRegistry, TrackHistogramSpec);
  };

  void process(FilteredCollision const& col, Tracks const& /*tracks*/)
  {
    ColHists.fill<femtounitedcolhistmanager::Mode::kANALYSIS>(col);
    auto TrackSlice = TrackPartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);

    for (auto const& track : TrackSlice) {
      TrackHists.fill<femtounitedtrackhistmanager::Mode::kANALYSIS>(track);
    };
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtounitedtrackqa>(cfgc),
  };
  return workflow;
}
