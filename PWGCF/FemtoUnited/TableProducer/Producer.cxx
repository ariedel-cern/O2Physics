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

#include <vector>
#include <string>

#include "CCDB/BasicCCDBManager.h"
#include "fairlogger/Logger.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/Expressions.h"
#include "Framework/Configurable.h"

// #include "EventFiltering/Zorro.h"

#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoVzerosDerived.h"

#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackTPCSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackTOFSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackTPCTOFSelection.h"
#include "PWGCF/FemtoUnited/Core/VzeroSelection.h"

#include "PWGCF/FemtoUnited/Utils/FemtoUtils.h"

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

namespace o2::analysis::femtounited
{
namespace consumedData
{
using Run3PpCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;
using Run3PpWithoutCentCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;

using Run3PpTracks =
  soa::Join<FullTracks, TracksDCA,
            pidTPCEl, pidTPCPi, pidTPCKa, pidTPCPr, pidTPCDe, pidTPCTr, pidTPCHe,
            pidTOFEl, pidTOFPi, pidTOFKa, pidTOFPr, pidTOFDe, pidTOFTr, pidTOFHe>;

using Run3PpQaFullPidTracks =
  soa::Join<FullTracks, TracksDCA,
            pidTPCFullEl, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTPCFullDe, pidTPCFullTr, pidTPCFullHe,
            pidTOFFullEl, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr, pidTOFFullDe, pidTOFFullTr, pidTOFFullHe,
            pidTOFbeta>;

using Run3PpVzeros = V0Datas;

} // namespace consumedData
} // namespace o2::analysis::femtounited

struct femtounitedProducer {

  // produced objectes
  Produces<o2::aod::FUCols> ProducedCollision;
  Produces<o2::aod::FUTracks> ProducedTracks;
  Produces<o2::aod::FUTrackMasks> ProducedTrackMasks;
  Produces<o2::aod::FUTrackDCAs> ProducedTrackDCAs;
  Produces<o2::aod::FUTrackExtras> ProducedTrackExtras;
  Produces<o2::aod::FUTrackPids> ProducedTrackPids;

  // configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("General");
    Configurable<std::string> CcdbUrl{"CccdbUrl", "http://alice-ccdb.cern.ch", "URL to ccdb"};
    Configurable<std::string> GrpPath{"GrpPath", "GLO/Config/GRPMagField", "Path to GRP object (Run3 -> GLO/Config/GRPMagField/Run2 -> GLO/GRP/GRP"};
  } ConfOptions;

  // Event selections
  struct : ConfigurableGroup {
    std::string prefix = std::string("CollisionSel");
    Configurable<float> ZvtxAbsMax{"ZvtxMax", 10.f, "Max. |z-Vertex| (cm)"};
    Configurable<bool> UseEventSel{"UseEventSel", true, "Use event selections (Run3 -> Sel8/Run2 -> Sel7)"};
  } ConfCollisionFilter;
  Filter CollisionFilter = nabs(o2::aod::collision::posZ) <= ConfCollisionFilter.ZvtxAbsMax;
  collisionselection::CollisionSelection CollisionSel;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackFilters");
    Configurable<float> PtMin{"PtMin", 0.15, "Minimum pT"};
    Configurable<float> PtMax{"PtMax", 6, "Maximum pT"};
    Configurable<float> EtaMin{"EtaMin", -0.8, "Minimum eta"};
    Configurable<float> EtaMax{"EtaMax", 0.8, "Maximum eta"};
    Configurable<float> PhiMin{"PhiMin", 0, "Minimum phi"};
    Configurable<float> PhiMax{"PhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  } ConfTrackFilters;

  Filter TrackFilter = track::pt >= ConfTrackFilters.PtMin && track::pt <= ConfTrackFilters.PtMax &&
                       track::eta >= ConfTrackFilters.EtaMin && track::eta <= ConfTrackFilters.EtaMax &&
                       track::phi >= ConfTrackFilters.PhiMin && track::phi <= ConfTrackFilters.PhiMax;

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
  trackselection::TrackSelection TrackSel;

  // track tpc bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackTPCBits");
    Configurable<std::vector<float>> Electron{"Electron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> Pion{"Pion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> Kaon{"Kaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> Proton{"Proton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> Deuteron{"Deuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> Triton{"Triton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> Helium{"Helium", {}, "Maximum |nsigma| for helium PID"};
  } ConfTrackTPCBits;
  tracktpcselection::TrackTpcSelection TrackTpcSel;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackTOFBits");
    Configurable<std::vector<float>> Electron{"Electron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> Pion{"Pion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> Kaon{"Kaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> Proton{"Proton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> Deuteron{"Deuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> Triton{"Triton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> Helium{"Helium", {}, "Maximum |nsigma| for helium3 PID"};
  } ConfTrackTOFBits;
  tracktofselection::TrackTofSelection TrackTofSel;

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
  tracktpctofselection::TrackTpcTofSelection TrackTpcTofSel;

  // v0 bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0Filters");
    Configurable<float> PtMin{"PtMin", 0.15, "Minimum pT"};
    Configurable<float> PtMax{"PtMax", 6, "Maximum pT"};
    Configurable<float> EtaMin{"EtaMin", -0.8, "Minimum eta"};
    Configurable<float> EtaMax{"EtaMax", 0.8, "Maximum eta"};
    Configurable<float> PhiMin{"PhiMin", 0, "Minimum phi"};
    Configurable<float> PhiMax{"PhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  } ConfVzeroFilters;
  Filter VzeroFilter = v0data::pt >= ConfVzeroFilters.PtMin && v0data::pt <= ConfVzeroFilters.PtMax &&
                       v0data::eta >= ConfVzeroFilters.EtaMin && v0data::eta <= ConfVzeroFilters.EtaMax &&
                       v0data::phi >= ConfVzeroFilters.PhiMin && v0data::phi <= ConfVzeroFilters.PhiMax;

  // v0 bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0Bits");
    Configurable<std::vector<float>> DcaDaughMax{"DcaDaughMax", {1.5}, "Maximum DCA between the daughters at decay vertex (cm)"};
    Configurable<std::vector<float>> CpaMin{"CpaMin", {0.99}, "Minimum cosine of pointing angle"};
    Configurable<std::vector<float>> TransRadMin{"TransRadMin", {0.2}, "Minimum transverse radius (cm)"};
    Configurable<std::vector<float>> TransRadMax{"TransRadMax", {100}, "Maximum transverse radius (cm)"};
    Configurable<std::vector<float>> DecayVtxMax{"DecayVtxMax", {100}, "Maximum distance in x,y,z of the decay vertex from primary vertex (cm)"};
  } ConfVzeroBits;
  vzeroselection::VzeroSelection VzeroSel;

  struct : ConfigurableGroup {
    std::string prefix = std::string("V0DaughterBits");
    Configurable<std::vector<float>> DcaMin{"DcaMin", {0.05}, "Minimum DCA of the daughters at primary vertex (cm)"};
    Configurable<std::vector<float>> TpcClustersMin{"TpcClustersMin", {70}, "Minimum number of TPC clusters for daughter tracks"};
    Configurable<std::vector<float>> TpcNsigmaMax{"TpcNsigmaMax", {5}, "Maximum |nsimga| TPC for daughter tracks"};
  } ConfVzeroDaughterBits;
  // vzeroselection::VzeroDaughterSelection VzeroDaughterSel;

  // histogramming
  HistogramRegistry hRegistry{"Producer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // data members
  int RunNumber = -1;
  float MagField = 0.f;
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB

  // functions
  void initFromCcdb(o2::aod::BCsWithTimestamps::iterator bc)
  {
    if (RunNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();

    static o2::parameters::GRPMagField* grpo = nullptr;
    grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ConfOptions.GrpPath.value, timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
      return;
    }
    MagField = 0.1 * grpo->getNominalL3Field(); // get magnetic field in tesla
    RunNumber = bc.runNumber();
  };

  template <modes::System sys, typename T1, typename T2>
  void fillCollision(T1 const& col, T2 const& tracks)
  {
    if constexpr (!modes::isSystemSet(sys, modes::System::kNoCentCal)) {
      ProducedCollision(col.posZ(),
                        col.multNTracksPV(),
                        col.centFT0M(),
                        utils::sphericity(tracks),
                        MagField);
    }

    if constexpr (modes::isSystemSet(sys, modes::System::kNoCentCal)) {
      ProducedCollision(col.posZ(),
                        col.multNTracksPV(),
                        0,
                        utils::sphericity(tracks),
                        MagField);
    }
  }

  template <modes::Mode mode, typename T>
  void fillTracks(T const& tracks)
  {
    for (const auto& track : tracks) {
      TrackSel.ApplySelections(track);
      if (TrackSel.getMinimalSelection()) {
        TrackTpcSel.ApplySelections(track);
        TrackTofSel.ApplySelections(track);
        TrackTpcTofSel.ApplySelections(track);
        if (TrackTpcSel.getAnySelection() || TrackTofSel.getAnySelection()) {
          if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
            ProducedTracks(ProducedCollision.lastIndex(), track.pt(), track.eta(), track.phi());
            ProducedTrackMasks(TrackSel.getBitmask(),
                               TrackTpcSel.getBitmask(),
                               TrackTofSel.getBitmask(),
                               TrackTpcTofSel.getBitmask());
          }

          if constexpr (modes::isModeSet(mode, modes::Mode::kQA)) {
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
  }

  template <modes::Mode mode, typename T>
  void fillV0s(T const& v0s)
  {
    for (const auto& v0 : v0s) {
      TrackSel.ApplySelections(v0);
      if (TrackSel.getMinimalSelection()) {
        TrackTpcSel.ApplySelections(v0);
        TrackTofSel.ApplySelections(v0);
        TrackTpcTofSel.ApplySelections(v0);
        if (TrackTpcSel.getAnySelection() || TrackTofSel.getAnySelection()) {
          if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
            ProducedTracks(ProducedCollision.lastIndex(), v0.pt(), v0.eta(), v0.phi());
            ProducedTrackMasks(TrackSel.getBitmask(),
                               TrackTpcSel.getBitmask(),
                               TrackTofSel.getBitmask(),
                               TrackTpcTofSel.getBitmask());
          }
        }
      }
    }
  }

  void init(InitContext&)
  {
    /// init ccdb
    ccdb->setURL(ConfOptions.CcdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    /// init collision selection
    CollisionSel.init(ConfCollisionFilter.UseEventSel.value);

    /// init track selections
    TrackSel.setCheckMinimalSelection(true);
    TrackSel.addSelection(ConfTrackBits.Sign.value, trackselection::kSign, limits::kEqual, false);
    TrackSel.addSelection(ConfTrackBits.TpcClustersMin.value, trackselection::kTPCnClsMin, limits::kLowerLimit, true);
    TrackSel.addSelection(ConfTrackBits.TpcCrossedRowsMin.value, trackselection::kTPCcRowsMin, limits::kLowerLimit, true);
    TrackSel.addSelection(ConfTrackBits.TpcSharedClustersMax.value, trackselection::kTPCsClsMax, limits::kUpperLimit, true);
    TrackSel.addSelection(ConfTrackBits.ItsClustersMin.value, trackselection::kITSnClsMin, limits::kLowerLimit, true);
    TrackSel.addSelection(ConfTrackBits.ItsIbClustersMin.value, trackselection::kITSnClsIbMin, limits::kLowerLimit, true);

    /// init track tpc selections
    TrackTpcSel.setCheckMinimalSelection(false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Electron.value, tracktpcselection::kElectron, limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Pion.value, tracktpcselection::kPion, limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Kaon.value, tracktpcselection::kKaon, limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Proton.value, tracktpcselection::kProton, limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Deuteron.value, tracktpcselection::kDeuteron, limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Triton.value, tracktpcselection::kTriton, limits::kAbsUpperLimit, false);
    TrackTpcSel.addSelection(ConfTrackTPCBits.Helium.value, tracktpcselection::kHelium, limits::kAbsUpperLimit, false);

    /// init track tof selections
    TrackTofSel.setCheckMinimalSelection(false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Electron.value, tracktofselection::kElectron, limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Pion.value, tracktofselection::kPion, limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Kaon.value, tracktofselection::kKaon, limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Proton.value, tracktofselection::kProton, limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Deuteron.value, tracktofselection::kDeuteron, limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Triton.value, tracktofselection::kTriton, limits::kAbsUpperLimit, false);
    TrackTofSel.addSelection(ConfTrackTOFBits.Helium.value, tracktofselection::kHelium, limits::kAbsUpperLimit, false);

    /// init track tpc+tof selections
    TrackTpcTofSel.setCheckMinimalSelection(false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Electron.value, tracktpctofselection::kElectron, limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Pion.value, tracktpctofselection::kPion, limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Kaon.value, tracktpctofselection::kKaon, limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Proton.value, tracktpctofselection::kProton, limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Deuteron.value, tracktpctofselection::kDeuteron, limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Triton.value, tracktpctofselection::kTriton, limits::kAbsUpperLimit, false);
    TrackTpcTofSel.addSelection(ConfTrackTPCTOFBits.Helium.value, tracktpctofselection::kHelium, limits::kAbsUpperLimit, false);
  }

  void processTracksRun3pp(Filtered<consumedData::Run3PpCollisions>::iterator const& col,
                           BCsWithTimestamps const&,
                           Filtered<consumedData::Run3PpTracks> const& tracks)
  {
    if (!CollisionSel.isSelected<modes::System::kPP_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());
    fillCollision<modes::System::kPP_Run3>(col, tracks);
    fillTracks<modes::Mode::kANALYSIS>(tracks);
  }
  PROCESS_SWITCH(femtounitedProducer, processTracksRun3pp, "Provide tracks for Run3 analysis", true);

  // void processVzerosRun3pp(Filtered<consumedData::Run3PpCollisions>::iterator const& col,
  //                          BCsWithTimestamps const&,
  //                          Filtered<consumedData::Run3PpTracks> const& tracks,
  //                          Filtered<consumedData::Run3PpVzeros> const& v0s)
  // {
  //   if (!CollisionSel.isSelected<modes::System::kPP_Run3>(col)) {
  //     return;
  //   }
  //   initFromCcdb(col.bc_as<BCsWithTimestamps>());
  //   fillCollision<modes::System::kPP_Run3>(col, tracks);
  //   fillV0s<modes::Mode::kANALYSIS>(v0s);
  // }
  // PROCESS_SWITCH(femtounitedProducer, processVzerosRun3pp, "Provide V0s for Run3 analysis", true);

  void proccessTracksRun3ppNoCent(Filtered<consumedData::Run3PpWithoutCentCollisions>::iterator const& col,
                                  BCsWithTimestamps const&,
                                  Filtered<consumedData::Run3PpTracks> const& tracks)
  {
    if (!CollisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());
    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    fillTracks<modes::Mode::kANALYSIS>(tracks);
  }
  PROCESS_SWITCH(femtounitedProducer, proccessTracksRun3ppNoCent, "Provide tracks for Run3 analysis (use when no centrality calibration is available)", false);

  void proccessQATracksRun3ppNoCent(Filtered<consumedData::Run3PpWithoutCentCollisions>::iterator const& col,
                                    BCsWithTimestamps const&,
                                    Filtered<consumedData::Run3PpQaFullPidTracks> const& tracks)
  {
    if (!CollisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());
    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    fillTracks<modes::Mode::kANALYSIS_QA>(tracks);
  }
  PROCESS_SWITCH(femtounitedProducer, proccessQATracksRun3ppNoCent, "Provide tracks with QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtounitedProducer>(cfgc)};
  return workflow;
}
