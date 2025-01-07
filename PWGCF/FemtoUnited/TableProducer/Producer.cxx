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
#include "Common/DataModel/PIDResponseITS.h"
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
// #include "PWGCF/FemtoUnited/Core/TrackTPCSelection.h"
// #include "PWGCF/FemtoUnited/Core/TrackTOFSelection.h"
// #include "PWGCF/FemtoUnited/Core/TrackTPCTOFSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackPidSelection.h"
#include "PWGCF/FemtoUnited/Core/VzeroSelection.h"
// #include "PWGCF/FemtoUnited/Core/VzeroDaughterSelection.h"

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

using Run3Tracks = soa::Join<Tracks, TracksExtra, TracksDCA,
                             pidTPCEl, pidTPCPi, pidTPCKa, pidTPCPr, pidTPCDe, pidTPCTr, pidTPCHe,
                             pidTOFEl, pidTOFPi, pidTOFKa, pidTOFPr, pidTOFDe, pidTOFTr, pidTOFHe>;

using Run3TracksFullPid =
  soa::Join<Tracks, TracksExtra, TracksDCA,
            pidTPCFullEl, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTPCFullDe, pidTPCFullTr, pidTPCFullHe,
            pidTOFFullEl, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr, pidTOFFullDe, pidTOFFullTr, pidTOFFullHe,
            pidTOFbeta>;

using Run3PpVzeros = V0Datas;

} // namespace consumedData
} // namespace o2::analysis::femtounited

struct femtounitedProducer {

  // produced objectes
  Produces<FUCols> ProducedCollision;

  Produces<FUTracks> ProducedTracks;
  Produces<FUTrackMasks> ProducedTrackMasks;
  Produces<FUTrackDCAs> ProducedTrackDCAs;
  Produces<FUTrackExtras> ProducedTrackExtras;
  Produces<FUTrackPids> ProducedTrackPids;

  Produces<FUVzeros> ProducedVzeros;
  Produces<FUVzeroMasks> ProducedVzeroMasks;
  Produces<FUVzeroExtras> ProducedVzeroExtras;
  Produces<FUVzeroDaus> ProducedVzeroDaus;
  Produces<FUVzeroDauExts> ProducedVzeroDauExts;

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
    Configurable<std::vector<string>> DcaxyMax{"DcaxyMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_xy| as a function of pT"};
    Configurable<std::vector<string>> DcazMax{"DcazMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_z| as a function of pT"};
  } ConfTrackBits;
  trackselection::TrackSelection TrackSel;

  // track pid bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackPidBits");
    Configurable<std::vector<float>> ItsElectron{"ItsElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> ItsPion{"ItsPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> ItsKaon{"ItsKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> ItsProton{"ItsProton", {}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> ItsDeuteron{"ItsDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> ItsTriton{"ItsTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> ItsHelium{"ItsHelium", {}, "Maximum |nsigma| for helium PID"};

    Configurable<std::vector<float>> TpcElectron{"TpcElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> TpcPion{"TpcPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> TpcKaon{"TpcKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> TpcProton{"TpcProton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> TpcDeuteron{"TpcDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> TpcTriton{"TpcTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> TpcHelium{"TpcHelium", {}, "Maximum |nsigma| for helium PID"};

    Configurable<std::vector<float>> TofElectron{"TofElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> TofPion{"TofPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> TofKaon{"TofKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> TofProton{"TofProton", {}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> TofDeuteron{"TofDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> TofTriton{"TofTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> TofHelium{"TofHelium", {}, "Maximum |nsigma| for helium PID"};

    Configurable<std::vector<float>> TpctofElectron{"TpctofElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> TpctofPion{"TpctofPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> TpctofKaon{"TpctofKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> TpctofProton{"TpctofProton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> TpctofDeuteron{"TpctofDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> TpctofTriton{"TpctofTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> TpctofHelium{"TpctofHelium", {}, "Maximum |nsigma| for helium PID"};

  } ConfTrackPidBits;
  trackpidselection::TrackPidSelection TrackPidSel;

  // v0 filters
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0Filters");
    Configurable<float> PtMin{"PtMin", 0, "Minimum pT"};
    Configurable<float> PtMax{"PtMax", 99, "Maximum pT"};
    Configurable<float> EtaMin{"EtaMin", -10, "Minimum eta"};
    Configurable<float> EtaMax{"EtaMax", 10, "Maximum eta"};
    Configurable<float> PhiMin{"PhiMin", 0, "Minimum phi"};
    Configurable<float> PhiMax{"PhiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  } ConfVzeroFilters;
  Filter VzeroFilter = v0data::pt >= ConfVzeroFilters.PtMin && v0data::pt <= ConfVzeroFilters.PtMax &&
                       v0data::eta >= ConfVzeroFilters.EtaMin && v0data::eta <= ConfVzeroFilters.EtaMax &&
                       v0data::phi >= ConfVzeroFilters.PhiMin && v0data::phi <= ConfVzeroFilters.PhiMax;

  // v0 bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0Bits");
    Configurable<std::vector<float>> Sign{"Sign", {-1, 1}, "+1 for particle and -1 for antiparticle"};
    Configurable<std::vector<float>> DcaDaughMax{"DcaDaughMax", {1.5}, "Maximum DCA between the daughters at decay vertex (cm)"};
    Configurable<std::vector<float>> CpaMin{"CpaMin", {0.99}, "Minimum cosine of pointing angle"};
    Configurable<std::vector<float>> TransRadMin{"TransRadMin", {0.2}, "Minimum transverse radius (cm)"};
    Configurable<std::vector<float>> TransRadMax{"TransRadMax", {100}, "Maximum transverse radius (cm)"};
    Configurable<std::vector<float>> DecayVtxMax{"DecayVtxMax", {100}, "Maximum distance in x,y,z of the decay vertex from primary vertex (cm)"};
    Configurable<float> KaonMassRejecionLow{"KaonMassRejecionLow", 0.48, "Lower limit of kaon mass (GeV/c^2) rejection (set negative to deactivate)"};
    Configurable<float> KaonMassRejecionHigh{"KaonMassRejecionUp", 0.515, "Upper limit of kaon mass (GeV/c^2) rejection (set negative to deactivate)"};
  } ConfVzeroBits;
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0DaughterBits");
    Configurable<std::vector<float>> DcaMin{"DcaMin", {0.05}, "Minimum DCA of the daughters from primary vertex (cm)"};
    Configurable<std::vector<float>> TpcClustersMin{"TpcClustersMin", {70}, "Minimum number of TPC clusters for daughter tracks"};
    Configurable<std::vector<float>> TpcNsigmaMax{"TpcNsigmaMax", {5}, "Maximum |nsimga| TPC for daughter tracks"};
  } ConfVzeroDaughterBits;
  vzeroselection::VzeroSelection VzeroSel;

  // histogramming
  HistogramRegistry hRegistry{"Producer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // data members
  int RunNumber = -1;
  float MagField = 0.f;
  Service<o2::ccdb::BasicCCDBManager> ccdb;          /// Accessing the CCDB
  std::vector<std::pair<int64_t, int64_t>> indexMap; // for mapping tracks to vzeros and more

  // functions
  void initFromCcdb(o2::aod::BCsWithTimestamps::iterator const bc)
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
  void fillTracks(T const& tracks, std::vector<std::pair<int64_t, int64_t>>& map)
  {
    for (const auto& track : tracks) {
      TrackSel.ApplySelections(track);
      if (TrackSel.getMinimalSelection()) {
        TrackPidSel.ApplySelections(track);
        if (TrackPidSel.getAnySelection()) {
          if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
            ProducedTracks(ProducedCollision.lastIndex(), track.pt(), track.eta(), track.phi());
            ProducedTrackMasks(TrackSel.getBitmask(), TrackPidSel.getBitmask());
          }

          if constexpr (modes::isModeSet(mode, modes::Mode::kQA)) {
            ProducedTrackDCAs(track.dcaXY(), track.dcaZ());
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
              track.itsNSigmaEl(),
              track.itsNSigmaPi(),
              track.itsNSigmaKa(),
              track.itsNSigmaPr(),
              track.itsNSigmaDe(),
              track.itsNSigmaTr(),
              track.itsNSigmaHe(),
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
          map.emplace_back(track.globalIndex(), ProducedTracks.lastIndex());
        }
      }
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fillV0s(T1 const& v0s, T2 const& tracks, std::vector<std::pair<int64_t, int64_t>> map)
  {
    for (const auto& v0 : v0s) {
      VzeroSel.ApplySelections(v0, tracks);
      if (VzeroSel.checkKaonMassLimit(v0) && VzeroSel.getMinimalSelection()) {
        auto posDaughter = v0.template posTrack_as<T2>();
        auto negDaughter = v0.template negTrack_as<T2>();
        if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
          ProducedVzeros(ProducedCollision.lastIndex(),
                         v0.pt(),
                         v0.eta(),
                         v0.phi(),
                         VzeroSel.getMass());
          ProducedVzeroMasks(VzeroSel.getBitmask());
          ProducedVzeroDaus(utils::getDaughterIndex(posDaughter.globalIndex(), map), posDaughter.pt(), posDaughter.eta(), posDaughter.phi(),
                            utils::getDaughterIndex(negDaughter.globalIndex(), map), negDaughter.pt(), negDaughter.eta(), negDaughter.phi());
        }
        if constexpr (modes::isModeSet(mode, modes::Mode::kQA)) {
          ProducedVzeroExtras(
            VzeroSel.getSign(),
            v0.dcaV0daughters(),
            v0.x(),
            v0.y(),
            v0.z(),
            v0.v0radius(),
            v0.mK0Short());
          ProducedVzeroDauExts(posDaughter.tpcNClsFound(), posDaughter.dcaXY(), posDaughter.dcaZ(), VzeroSel.getPosDaughterTpcNsigma(),
                               negDaughter.tpcNClsFound(), negDaughter.dcaXY(), negDaughter.dcaZ(), VzeroSel.getNegDaughterTpcNsigma());
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
    TrackSel.addSelection(ConfTrackBits.DcaxyMax.name, ConfTrackFilters.PtMin.value, ConfTrackFilters.PtMax.value, ConfTrackBits.DcaxyMax.value, trackselection::kDCAxyMax, limits::kAbsUpperLimit, true);
    TrackSel.addSelection(ConfTrackBits.DcazMax.name, ConfTrackFilters.PtMin.value, ConfTrackFilters.PtMax.value, ConfTrackBits.DcazMax.value, trackselection::kDCAzMax, limits::kAbsUpperLimit, true);

    /// init track pid selections
    TrackPidSel.setCheckMinimalSelection(false);

    TrackPidSel.addSelection(ConfTrackPidBits.ItsElectron.value, trackpidselection::kItsElectron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.ItsPion.value, trackpidselection::kItsPion, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.ItsKaon.value, trackpidselection::kItsKaon, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.ItsProton.value, trackpidselection::kItsProton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.ItsDeuteron.value, trackpidselection::kItsDeuteron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.ItsTriton.value, trackpidselection::kItsTriton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.ItsHelium.value, trackpidselection::kItsHelium, limits::kAbsUpperLimit, false);

    TrackPidSel.addSelection(ConfTrackPidBits.TpcElectron.value, trackpidselection::kTpcElectron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpcPion.value, trackpidselection::kTpcPion, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpcKaon.value, trackpidselection::kTpcKaon, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpcProton.value, trackpidselection::kTpcProton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpcDeuteron.value, trackpidselection::kTpcDeuteron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpcTriton.value, trackpidselection::kTpcTriton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpcHelium.value, trackpidselection::kTpcHelium, limits::kAbsUpperLimit, false);

    TrackPidSel.addSelection(ConfTrackPidBits.TofElectron.value, trackpidselection::kTofElectron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TofPion.value, trackpidselection::kTofPion, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TofKaon.value, trackpidselection::kTofKaon, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TofProton.value, trackpidselection::kTofProton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TofDeuteron.value, trackpidselection::kTofDeuteron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TofTriton.value, trackpidselection::kTofTriton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TofHelium.value, trackpidselection::kTofHelium, limits::kAbsUpperLimit, false);

    TrackPidSel.addSelection(ConfTrackPidBits.TpctofElectron.value, trackpidselection::kTpctofElectron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpctofPion.value, trackpidselection::kTpctofPion, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpctofKaon.value, trackpidselection::kTpctofKaon, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpctofProton.value, trackpidselection::kTpctofProton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpctofDeuteron.value, trackpidselection::kTpctofDeuteron, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpctofTriton.value, trackpidselection::kTpctofTriton, limits::kAbsUpperLimit, false);
    TrackPidSel.addSelection(ConfTrackPidBits.TpctofHelium.value, trackpidselection::kTpctofHelium, limits::kAbsUpperLimit, false);

    /// init vzero selections
    VzeroSel.setCheckMinimalSelection(true);
    VzeroSel.addSelection(ConfVzeroBits.Sign.value, vzeroselection::kSign, limits::kEqual, false);
    VzeroSel.addSelection(ConfVzeroBits.DcaDaughMax.value, vzeroselection::kDcaDaughMax, limits::kAbsUpperLimit, true);
    VzeroSel.addSelection(ConfVzeroBits.CpaMin.value, vzeroselection::kCpaMin, limits::kLowerLimit, true);
    VzeroSel.addSelection(ConfVzeroBits.TransRadMin.value, vzeroselection::kTransRadMin, limits::kLowerLimit, true);
    VzeroSel.addSelection(ConfVzeroBits.TransRadMax.value, vzeroselection::kTransRadMax, limits::kUpperLimit, true);
    // vzero positiv daughter selections
    VzeroSel.addSelection(ConfVzeroDaughterBits.DcaMin.value, vzeroselection::kPosDauDcaMin, limits::kLowerLimit, true);
    VzeroSel.addSelection(ConfVzeroDaughterBits.TpcClustersMin.value, vzeroselection::kPosDauTpcClsMin, limits::kLowerLimit, true);
    VzeroSel.addSelection(ConfVzeroDaughterBits.TpcNsigmaMax.value, vzeroselection::kPosDauTpcNsigmaMax, limits::kAbsUpperLimit, true);
    // vzero negative daughter selections
    VzeroSel.setCheckMinimalSelection(true);
    VzeroSel.addSelection(ConfVzeroDaughterBits.DcaMin.value, vzeroselection::kNegDauDcaMin, limits::kLowerLimit, true);
    VzeroSel.addSelection(ConfVzeroDaughterBits.TpcClustersMin.value, vzeroselection::kNegDauTpcClsMin, limits::kLowerLimit, true);
    VzeroSel.addSelection(ConfVzeroDaughterBits.TpcNsigmaMax.value, vzeroselection::kNegDauTpcNsigmaMax, limits::kAbsUpperLimit, true);
    // kaon mass limits
    VzeroSel.setKaonMassLimits(ConfVzeroBits.KaonMassRejecionLow.value, ConfVzeroBits.KaonMassRejecionHigh.value);
  }

  // proccess functions
  // produce tracks for analysis
  void processTracksRun3pp(Filtered<consumedData::Run3PpCollisions>::iterator const& col,
                           BCsWithTimestamps const&,
                           Filtered<consumedData::Run3Tracks> const& tracks)
  {
    if (!CollisionSel.isSelected<modes::System::kPP_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumedData::Run3Tracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS>(tracksWithItsPid, indexMap);
  }
  PROCESS_SWITCH(femtounitedProducer, processTracksRun3pp, "Provide tracks for Run3 analysis", true);

  // produce tracks for analysis (without centrality)
  void proccessTracksRun3ppNoCent(Filtered<consumedData::Run3PpWithoutCentCollisions>::iterator const& col,
                                  BCsWithTimestamps const&,
                                  Filtered<consumedData::Run3Tracks> const& tracks)
  {
    if (!CollisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumedData::Run3Tracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS>(tracksWithItsPid, indexMap);
  }
  PROCESS_SWITCH(femtounitedProducer, proccessTracksRun3ppNoCent, "Provide tracks for Run3 analysis (use when no centrality calibration is available)", false);

  // produce tracks for QA (without centrality)
  void proccessQaTracksRun3ppNoCent(Filtered<consumedData::Run3PpWithoutCentCollisions>::iterator const& col,
                                    BCsWithTimestamps const&,
                                    Filtered<consumedData::Run3TracksFullPid> const& tracks)
  {
    if (!CollisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumedData::Run3TracksFullPid, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS_QA>(tracksWithItsPid, indexMap);
  }
  PROCESS_SWITCH(femtounitedProducer, proccessQaTracksRun3ppNoCent, "Provide tracks for Run2 with QA", false);

  // produce tracks and v0s for QA (without centrality)
  void processQaTracksVzerosRun3ppNoCent(Filtered<consumedData::Run3PpWithoutCentCollisions>::iterator const& col,
                                         BCsWithTimestamps const&,
                                         Filtered<consumedData::Run3TracksFullPid> const& tracks,
                                         Filtered<consumedData::Run3PpVzeros> const& v0s)
  {
    if (!CollisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = Attach<consumedData::Run3TracksFullPid, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                   pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS_QA>(tracksWithItsPid, indexMap);
    fillV0s<modes::Mode::kANALYSIS_QA>(v0s, tracks, indexMap);
  }
  PROCESS_SWITCH(femtounitedProducer, processQaTracksVzerosRun3ppNoCent, "Provide Tracks and V0s for Run3 with QA (no centrality calibration)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtounitedProducer>(cfgc)};
  return workflow;
}
