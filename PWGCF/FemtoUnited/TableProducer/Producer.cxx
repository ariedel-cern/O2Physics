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

/// \file femtoDreamProducerTask.cxx
/// \brief Tasks that produces the track tables used for the pairing
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <CCDB/BasicCCDBManager.h>
#include <fairlogger/Logger.h>
#include <vector>
#include <string>

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Expressions.h"
#include "Framework/Configurable.h"

// #include "EventFiltering/Zorro.h"

#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/Core/TrackSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackTPCSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackTOFSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackTPCTOFSelection.h"

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

namespace o2::aod::femtounited
{
namespace consumedData
{
using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;
using Tracks =
  soa::Join<aod::FullTracks, aod::TracksDCA,
            aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
            aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
            aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa,
            aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe,
            aod::pidTOFbeta, aod::pidTOFmass>;

} // namespace consumedData
} // namespace o2::aod::femtounited

struct femtounitedProducer {

  // produce objectes
  Produces<o2::aod::FUCols> ProducedCollision;
  Produces<o2::aod::FUTracks> ProducedTracks;
  Produces<o2::aod::FUTrackMasks> ProducedTrackMasks;
  Produces<o2::aod::FUTrackDCAs> ProducedTrackDCAs;
  Produces<o2::aod::FUTrackExtras> ProducedTrackExtras;
  Produces<o2::aod::FUTrackPids> ProducedTrackPids;

  // configurables
  // Event selections
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<float> ZvtxAbsMax{"ZvtxMax", 10.f, "Max. |z-Vertex| (cm)"};
    Configurable<bool> UseSel8{"UseSel8", true, "Use sel8 selections"};
  } ConfEventFilters;

  Filter EventFilter = nabs(o2::aod::collision::posZ) <= ConfEventFilters.ZvtxAbsMax;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackFilters");
    Configurable<float> PtMin{"PtMin", 0.15, "Minimum pT"};
    Configurable<float> PtMax{"PtMax", 6, "Maximum pT"};
    Configurable<float> EtaMin{"EtaMin", -0.8, "Minimum eta"};
    Configurable<float> EtaMax{"EtaMax", 0.8, "Maximum eta"};
    Configurable<float> PhiMin{"PhiMin", 0, "Minimum phi"};
    Configurable<float> PhiMax{"PhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  } TrackFilters;

  Filter TrackFilter = track::pt >= TrackFilters.PtMin && track::pt <= TrackFilters.PtMax &&
                       track::eta >= TrackFilters.EtaMin && track::eta <= TrackFilters.EtaMax &&
                       track::phi >= TrackFilters.PhiMin && track::phi <= TrackFilters.PhiMax;

  // track bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackBits");
    Configurable<std::vector<float>> Sign{"Sign", {1, -1}, "Sign (|Charge|) of the track"};
    Configurable<std::vector<float>> TpcClustersMin{"TpcClustersMin", {80}, "Minimum number of clusters in TPC"};
    Configurable<std::vector<float>> TpcCrossedRowsMin{"TpcCrossedRowsMin", {70}, "Minimum number of crossed rows in TPC"};
    Configurable<std::vector<float>> TpcSharedClustersMax{"TpcSharedClustersMax", {2}, "Maximum number of shared clusters in TPC"};
    Configurable<std::vector<float>> ItsClustersMin{"ItsClustersMin", {7}, "Minimum number of clusters in ITS"};
    Configurable<std::vector<float>> ItsIbClustersMin{"ItsIbClustersMin", {3}, "Minimum number of clusters in inner barrel (max 3) of ITS"};
  } ConfTrackBits;
  o2::analysis::femtounited::TrackSelection TrackSel;

  // track tpc bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackTPCBits");
    Configurable<std::vector<float>> Electron{"Electron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> Pion{"Pion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> Kaon{"Kaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> Proton{"Proton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> Deuteron{"Deuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> Triton{"Triton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> Helium3{"Helium3", {}, "Maximum |nsigma| for helium3 PID"};
  } ConfTrackTPCBits;
  o2::analysis::femtounited::TrackTpcSelection TrackTpcSel;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackTOFBits");
    Configurable<std::vector<float>> Electron{"Electron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> Pion{"Pion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> Kaon{"Kaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> Proton{"Proton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> Deuteron{"Deuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> Triton{"Triton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> Helium3{"Helium3", {}, "Maximum |nsigma| for helium3 PID"};
  } ConfTrackTOFBits;
  o2::analysis::femtounited::TrackTpcSelection TrackTofSel;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackTPCTOFBits");
    Configurable<std::vector<float>> Electron{"Electron", {}, "Maximum (nsigma_TPC^2+nsigma_TOF^2)^(1/2) for electron PID"};
    Configurable<std::vector<float>> Pion{"Pion", {}, "Maximum (nsigma_TPC^2+nsigma_TOF^2)^(1/2) for pion PID"};
    Configurable<std::vector<float>> Kaon{"Kaon", {}, "Maximum (nsigma_TPC^2+nsigma_TOF^2)^(1/2) for kaon PID"};
    Configurable<std::vector<float>> Proton{"Proton", {3}, "Maximum (nsigma_TPC^2+nsigma_TOF^2)^(1/2) for proton PID"};
    Configurable<std::vector<float>> Deuteron{"Deuteron", {}, "Maximum (nsigma_TPC^2+nsigma_TOF^2)^(1/2) for deuteron PID"};
    Configurable<std::vector<float>> Triton{"Triton", {}, "Maximum (nsigma_TPC^2+nsigma_TOF^2)^(1/2) for trition PID"};
    Configurable<std::vector<float>> Helium{"Helium", {}, "Maximum (nsigma_TPC^2+nsigma_TOF^2)^(1/2) for helium PID"};
  } ConfTrackTPCTOFBits;
  o2::analysis::femtounited::TrackTpcTofSelection TrackTpcTofSel;

  // histogramming
  HistogramRegistry hRegistry{"Producer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // ccdb
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  void init(InitContext&)
  {
    /// init ccdb
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    /// init track selections
    TrackSel.setCheckMinimalSelection(true);
    TrackSel.addSelection(ConfTrackBits.Sign.value, TrackSel::kSign, Limits::kEqual, false);
    TrackSel.addSelection(ConfTrackBits.TpcClustersMin.value, TrackSel::kTPCnClsMin, Limits::kLowerLimit, true);
    TrackSel.addSelection(ConfTrackBits.TpcCrossedRowsMin.value, TrackSel::kTPCcRowsMin, Limits::kLowerLimit, true);
    TrackSel.addSelection(ConfTrackBits.TpcSharedClustersMax.value, TrackSel::kTPCsClsMax, Limits::kUpperLimit, true);
    TrackSel.addSelection(ConfTrackBits.ItsClustersMin.value, TrackSel::kITSnClsMin, Limits::kLowerLimit, true);
    TrackSel.addSelection(ConfTrackBits.ItsIbClustersMin.value, TrackSel::kITSnClsIbMin, Limits::kLowerLimit, true);

    /// init track tpc selections
    TrackTpcSel.setCheckMinimalSelection(false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Electron.value, TrackTpcSel::kElectron, Limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Pion.value, TrackTpcSel::kPion, Limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Kaon.value, TrackTpcSel::kKaon, Limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Proton.value, TrackTpcSel::kProton, Limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Deuteron.value, TrackTpcSel::kDeuteron, Limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Triton.value, TrackTpcSel::kTriton, Limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Helium3.value, TrackTpcSel::kHelium3, Limits::kAbsUpperLimit, false);

    /// init track tof selections
    TrackTofSel.setCheckMinimalSelection(false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Electron.value, TrackTofSel::kElectron, Limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Pion.value, TrackTofSel::kPion, Limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Kaon.value, TrackTofSel::kKaon, Limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Proton.value, TrackTofSel::kProton, Limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Deuteron.value, TrackTofSel::kDeuteron, Limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Triton.value, TrackTofSel::kTriton, Limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Helium3.value, TrackTofSel::kHelium3, Limits::kAbsUpperLimit, false);

    /// init track tpc+tof selections
    TrackTpcTofSel.setCheckMinimalSelection(false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Electron.value, TrackTpcTofSel::kElectron, Limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Pion.value, TrackTpcTofSel::kPion, Limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Kaon.value, TrackTpcTofSel::kKaon, Limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Proton.value, TrackTpcTofSel::kProton, Limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Deuteron.value, TrackTpcTofSel::kDeuteron, Limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Triton.value, TrackTpcTofSel::kTriton, Limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Helium.value, TrackTpcTofSel::kHelium3, Limits::kAbsUpperLimit, false);
  }

  void processTracks(Filtered<femtounited::consumedData::Collisions>::iterator const& col,
                     o2::aod::BCsWithTimestamps const&,
                     Filtered<femtounited::consumedData::Tracks> const& tracks)
  {
    if (ConfEventFilters.UseSel8.value == true && col.sel8() == false) {
      return;
    }

    ProducedCollision(col.posZ(),
                      col.multNTracksPV(),
                      col.centFT0M(),
                      1);

    for (const auto& track : tracks) {
      TrackSel.ApplySelections(track);
      if (TrackSel.getMinimalSelection()) {
        TrackTpcSel.ApplySelections(track);
        TrackTofSel.ApplySelections(track);
        TrackTpcTofSel.ApplySelections(track);
        if (TrackTpcSel.getAnySelection() || TrackTpcTofSel.getAnySelection()) {
          ProducedTracks(ProducedCollision.lastIndex(), track.pt(), track.eta(), track.phi());
          ProducedTrackMasks(TrackSel.getBitmask(),
                             TrackTpcSel.getBitmask(),
                             TrackTofSel.getBitmask(),
                             TrackTpcTofSel.getBitmask());
        }
      }
    }
  }
  PROCESS_SWITCH(femtounitedProducer, processTracks, "Provide tracks", true);

  void processTrackQA(
    Filtered<femtounited::consumedData::Collisions>::iterator const& col,
    Filtered<femtounited::consumedData::Tracks> const& tracks)
  {
    if (ConfEventFilters.UseSel8.value == true && col.sel8() == false) {
      return;
    }
    ProducedCollision(col.posZ(), 1, 1, 1);
    for (const auto& track : tracks) {
      TrackSel.ApplySelections(track);
      if (TrackSel.getMinimalSelection()) {
        TrackTpcSel.ApplySelections(track);
        TrackTofSel.ApplySelections(track);
        TrackTpcTofSel.ApplySelections(track);
        if (TrackTpcSel.getAnySelection() || TrackTofSel.getAnySelection() || TrackTpcTofSel.getAnySelection()) {
          ProducedTracks(ProducedCollision.lastIndex(), track.pt(), track.eta(), track.phi());
          ProducedTrackMasks(TrackSel.getBitmask(),
                             TrackTpcSel.getBitmask(),
                             TrackTofSel.getBitmask(),
                             TrackTpcTofSel.getBitmask());
          ProducedTrackDCAs(track.dcaXY(),
                            track.dcaZ());
          ProducedTrackExtras(track.sign(),
                              track.isPVContributor(),
                              track.itsNCls(),
                              track.itsNClsInnerBarrel(),
                              track.itsChi2NCl(),
                              track.itsClusterSizes(),
                              track.tpcSignal(),
                              track.tpcInnerParam(),
                              track.tpcNClsFound(),
                              track.tpcNClsCrossedRows(),
                              track.tpcNClsShared(),
                              track.beta());
          ProducedTrackPids(
            track.tpcNSigmaEl(),
            track.tpcNSigmaPi(),
            track.tpcNSigmaKa(),
            track.tpcNSigmaPr(),
            track.tpcNSigmaDe(),
            track.tpcNSigmaTr(),
            track.tpcNSigmaHe(),
            track.tofNSigmaEl(),
            track.tofNSigmaPi(),
            track.tofNSigmaKa(),
            track.tofNSigmaPr(),
            track.tofNSigmaDe(),
            track.tofNSigmaTr(),
            track.tofNSigmaHe());
        }
      }
    }
  }
  PROCESS_SWITCH(femtounitedProducer, processTrackQA, "Provide tracks with QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtounitedProducer>(cfgc)};
  return workflow;
}
