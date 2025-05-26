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

/// \file FemtoUnitedProducer.cxx
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
#include "PWGCF/FemtoUnited/Core/TrackPidSelection.h"
#include "PWGCF/FemtoUnited/Core/VzeroSelection.h"
#include "PWGCF/FemtoUnited/Core/VzeroDaughterPidSelection.h"

#include "PWGCF/FemtoUnited/Utils/FemtoUtils.h"

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

namespace o2::analysis::femtounited
{
namespace consumeddata
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

} // namespace consumeddata
} // namespace o2::analysis::femtounited

struct FemtoUnitedProducer {

  // produced objectes
  Produces<FUCols> producedCollision;

  Produces<FUTracks> producedTracks;
  Produces<FUTrackMasks> producedTrackMasks;
  Produces<FUTrackDCAs> producedTrackDCAs;
  Produces<FUTrackExtras> producedTrackExtras;
  Produces<FUTrackPids> producedTrackPids;

  Produces<FUVzeros> producedVzeros;
  Produces<FUVzeroMasks> producedVzeroMasks;
  Produces<FUVzeroExtras> producedVzeroExtras;
  Produces<FUVzeroDaus> producedVzeroDaus;
  Produces<FUVzeroDauExts> producedVzeroDauExts;

  // configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("General");
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL to ccdb"};
    Configurable<std::string> grpPath{"grpPath", "GLO/Config/GRPMagField", "Path to GRP object (Run3 -> GLO/Config/GRPMagField/Run2 -> GLO/GRP/GRP"};
  } ConfOptions;

  // Event selections
  struct : ConfigurableGroup {
    std::string prefix = std::string("CollisionSel");
    Configurable<float> zVtxAbsMax{"zVtxAbsMax", 10.f, "Max. |z-Vertex| (cm)"};
    Configurable<bool> useEventSel{"useEventSel", true, "Use event selections (Run3 -> Sel8/Run2 -> Sel7)"};
  } ConfCollisionFilter;
  Filter collisionFilter = nabs(o2::aod::collision::posZ) <= ConfCollisionFilter.zVtxAbsMax;
  collisionselection::CollisionSelection collisionSel;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackFilters");
    Configurable<float> ptMin{"ptMin", 0.15, "Minimum pT"};
    Configurable<float> ptMax{"ptMax", 6, "Maximum pT"};
    Configurable<float> etaMin{"etaMin", -0.8, "Minimum eta"};
    Configurable<float> etaMax{"etaMax", 0.8, "Maximum eta"};
    Configurable<float> phiMin{"phiMin", 0, "Minimum phi"};
    Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  } ConfTrackFilters;

  Filter trackFilter = track::pt >= ConfTrackFilters.ptMin && track::pt <= ConfTrackFilters.ptMax &&
                       track::eta >= ConfTrackFilters.etaMin && track::eta <= ConfTrackFilters.etaMax &&
                       track::phi >= ConfTrackFilters.phiMin && track::phi <= ConfTrackFilters.phiMax;

  // track bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackBits");
    Configurable<bool> fillPositiveTracks{"fillPositiveTracks", true, "Fill if track if it has postivie charge"};
    Configurable<bool> fillNegativeTracks{"fillNegativeTracks", true, "Fill if track if it has negative charge"};
    Configurable<std::vector<float>> tpcClustersMin{"tpcClustersMin", {80}, "Minimum number of clusters in TPC"};
    Configurable<std::vector<float>> tpcCrossedRowsMin{"tpcCrossedRowsMin", {70}, "Minimum number of crossed rows in TPC"};
    Configurable<std::vector<float>> tpcSharedClustersMax{"tpcSharedClustersMax", {2}, "Maximum number of shared clusters in TPC"};
    Configurable<std::vector<float>> itsClustersMin{"itsClustersMin", {7}, "Minimum number of clusters in ITS"};
    Configurable<std::vector<float>> itsIbClustersMin{"itsIbClustersMin", {3}, "Minimum number of clusters in inner barrel (max 3) of ITS"};
    Configurable<std::vector<string>> dcaxyMax{"dcaxyMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_xy| as a function of pT"};
    Configurable<std::vector<string>> dcazMax{"dcazMax", {"0.0105+(0.035/x^(1.1))"}, "Maximum |dca_z| as a function of pT"};
  } ConfTrackBits;
  trackselection::TrackSelection trackSel;

  // track pid bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackPidBits");
    Configurable<std::vector<float>> itsElectron{"itsElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> itsPion{"itsPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> itsKaon{"itsKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> itsProton{"itsProton", {}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> itsDeuteron{"itsDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> itsTriton{"itsTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> itsHelium{"itsHelium", {}, "Maximum |nsigma| for helium PID"};

    Configurable<std::vector<float>> tpcElectron{"tpcElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> tpcPion{"tpcPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> tpcKaon{"tpcKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> tpcProton{"tpcProton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> tpcDeuteron{"tpcDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> tpcTriton{"tpcTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> tpcHelium{"tpcHelium", {}, "Maximum |nsigma| for helium PID"};

    Configurable<std::vector<float>> tofElectron{"tofElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> tofPion{"tofPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> tofKaon{"tofKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> tofProton{"tofProton", {}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> tofDeuteron{"tofDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> tofTriton{"tofTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> tofHelium{"tofHelium", {}, "Maximum |nsigma| for helium PID"};

    Configurable<std::vector<float>> tpctofElectron{"tpctofElectron", {}, "Maximum |nsigma| for electron PID"};
    Configurable<std::vector<float>> tpctofPion{"tpctofPion", {}, "Maximum |nsigma| for pion PID"};
    Configurable<std::vector<float>> tpctofKaon{"tpctofKaon", {}, "Maximum |nsigma| for kaon PID"};
    Configurable<std::vector<float>> tpctofProton{"tpctofProton", {3}, "Maximum |nsigma| for proton PID"};
    Configurable<std::vector<float>> tpctofDeuteron{"tpctofDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
    Configurable<std::vector<float>> tpctofTriton{"tpctofTriton", {}, "Maximum |nsigma| for trition PID"};
    Configurable<std::vector<float>> tpctofHelium{"tpctofHelium", {}, "Maximum |nsigma| for helium PID"};

  } ConfTrackPidBits;
  trackpidselection::TrackPidSelection trackPidSel;

  // v0 filters
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0Filters");
    Configurable<float> ptMin{"ptMin", 0, "Minimum pT"};
    Configurable<float> ptMax{"ptMax", 99, "Maximum pT"};
    Configurable<float> etaMin{"etaMin", -10, "Minimum eta"};
    Configurable<float> etaMax{"etaMax", 10, "Maximum eta"};
    Configurable<float> phiMin{"phiMin", 0, "Minimum phi"};
    Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  } ConfVzeroFilters;
  Filter vzeroFilter = v0data::pt >= ConfVzeroFilters.ptMin && v0data::pt <= ConfVzeroFilters.ptMax &&
                       v0data::eta >= ConfVzeroFilters.etaMin && v0data::eta <= ConfVzeroFilters.etaMax &&
                       v0data::phi >= ConfVzeroFilters.phiMin && v0data::phi <= ConfVzeroFilters.phiMax;

  // v0 bits
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0Bits");
    Configurable<std::vector<float>> dcaDaughMax{"dcaDaughMax", {1.5}, "Maximum DCA between the daughters at decay vertex (cm)"};
    Configurable<std::vector<float>> cpaMin{"cpaMin", {0.99}, "Minimum cosine of pointing angle"};
    Configurable<std::vector<float>> transRadMin{"transRadMin", {0.2}, "Minimum transverse radius (cm)"};
    Configurable<std::vector<float>> transRadMax{"transRadMax", {100}, "Maximum transverse radius (cm)"};
    Configurable<std::vector<float>> decayVtxMax{"decayVtxMax", {100}, "Maximum distance in x,y,z of the decay vertex from primary vertex (cm)"};
    Configurable<float> kaonMassRejectionLow{"kaonMassRejectionLow", 0.48, "Lower limit of kaon mass (GeV/c^2) rejection (set negative to deactivate)"};
    Configurable<float> kaonMassRejectionHigh{"kaonMassRejectionHigh", 0.515, "Upper limit of kaon mass (GeV/c^2) rejection (set negative to deactivate)"};
  } ConfVzeroBits;
  struct : ConfigurableGroup {
    std::string prefix = std::string("V0DaughterBits");
    Configurable<std::vector<float>> dcaMin{"dcaMin", {0.05}, "Minimum DCA of the daughters from primary vertex (cm)"};
    Configurable<std::vector<float>> tpcClustersMin{"tpcClustersMin", {70}, "Minimum number of TPC clusters for daughter tracks"};
    Configurable<std::vector<float>> posDaughProtonNsigmaMax{"posDaughProtonNsigmaMax", {5}, "Maximum |nsimga_Proton| TPC for positive daughter tracks"};
    Configurable<std::vector<float>> posDaughPionNsigmaMax{"posDaughPionNsigmaMax", {5}, "Maximum |nsimga_Pion| TPC for positive daughter tracks"};
    Configurable<std::vector<float>> negDaughProtonNsigmaMax{"negDaughProtonNsigmaMax", {5}, "Maximum |nsimga_Proton| TPC negative for daughter tracks"};
    Configurable<std::vector<float>> negDaughPionNsigmaMax{"negDaughPionNsigmaMax", {5}, "Maximum |nsimga_Pion| TPC for negative daughter tracks"};
  } ConfVzeroDaughterBits;
  vzeroselection::VzeroSelection vzeroSel;
  vzerodaughterpidselection::VzeroDaughterPidSelection vzeroDaugherPidSel;

  // histogramming
  HistogramRegistry hRegistry{"Producer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // data members
  int runNumber = -1;
  float magField = 0.f;
  Service<o2::ccdb::BasicCCDBManager> ccdb;      /// Accessing the CCDB
  std::unordered_map<int64_t, int64_t> indexMap; // for mapping tracks to vzeros and more

  // functions
  void initFromCcdb(o2::aod::BCsWithTimestamps::iterator const bc)
  {
    if (runNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();

    static o2::parameters::GRPMagField* grpo = nullptr;
    grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ConfOptions.grpPath.value, timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
      return;
    }
    magField = 0.1 * grpo->getNominalL3Field(); // get magnetic field in tesla
    runNumber = bc.runNumber();
  };

  template <modes::System sys, typename T1, typename T2>
  void fillCollision(T1 const& col, T2 const& tracks)
  {
    if constexpr (!modes::isSystemSet(sys, modes::System::kNoCentCal)) {
      producedCollision(col.posZ(),
                        col.multNTracksPV(),
                        col.centFT0M(),
                        utils::sphericity(tracks),
                        magField);
    }

    if constexpr (modes::isSystemSet(sys, modes::System::kNoCentCal)) {
      producedCollision(col.posZ(),
                        col.multNTracksPV(),
                        0,
                        utils::sphericity(tracks),
                        magField);
    }
  }

  template <modes::Mode mode, typename T, typename I>
  void fillTracks(T const& tracks, std::unordered_map<I, I>& map)
  {
    for (const auto& track : tracks) {
      trackSel.applySelections(track);
      if (trackSel.getMinimalSelection()) {
        trackPidSel.applySelections(track);
        if (trackPidSel.getAnySelection()) {
          if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
            bool hasPositiveCharge = track.sign() > 0;
            if ((ConfTrackBits.fillPositiveTracks && hasPositiveCharge) || (ConfTrackBits.fillNegativeTracks && !hasPositiveCharge)) {
              producedTracks(producedCollision.lastIndex(), track.pt(), track.eta(), track.phi(), hasPositiveCharge);
              producedTrackMasks(trackSel.getBitmask(), trackPidSel.getBitmask());
            }
          }

          if constexpr (modes::isModeSet(mode, modes::Mode::kQA)) {
            producedTrackDCAs(track.dcaXY(), track.dcaZ());
            producedTrackExtras(track.sign(),
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
            producedTrackPids(
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
          // map.emplace_back(track.globalIndex(), producedTracks.lastIndex());
          map[track.globalIndex()] = producedTracks.lastIndex();
        }
      }
    }
  }

  template <modes::Mode mode, typename V, typename T, typename I>
  void fillV0s(V const& v0s, T const& tracks, std::unordered_map<I, I>& map)
  {
    for (const auto& v0 : v0s) {
      vzeroSel.applySelections(v0, tracks);
      if (vzeroSel.checkKaonMassLimit(v0) && vzeroSel.getMinimalSelection()) {
        auto posDaughter = v0.template posTrack_as<T>();
        auto negDaughter = v0.template negTrack_as<T>();

        vzeroDaugherPidSel.applySelections(posDaughter, negDaughter);

        if (vzeroDaugherPidSel.isLambda() || vzeroDaugherPidSel.isAntiLambda()) {
          if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
            producedVzeros(producedCollision.lastIndex(),
                           v0.pt(),
                           v0.eta(),
                           v0.phi(),
                           v0.mLambda(),
                           v0.mAntiLambda());
            producedVzeroMasks(vzeroSel.getBitmask(), vzeroDaugherPidSel.getBitmask());
            producedVzeroDaus(utils::getDaughterIndex(posDaughter.globalIndex(), map), posDaughter.pt(), posDaughter.eta(), posDaughter.phi(),
                              utils::getDaughterIndex(negDaughter.globalIndex(), map), negDaughter.pt(), negDaughter.eta(), negDaughter.phi());
          }
          if constexpr (modes::isModeSet(mode, modes::Mode::kQA)) {
            producedVzeroExtras(
              v0.dcaV0daughters(),
              v0.x(),
              v0.y(),
              v0.z(),
              v0.v0radius(),
              v0.mK0Short());
            producedVzeroDauExts(posDaughter.tpcNClsFound(), posDaughter.dcaXY(), posDaughter.dcaZ(), posDaughter.tpcNSigmaPr(), posDaughter.tpcNSigmaPi(),
                                 negDaughter.tpcNClsFound(), negDaughter.dcaXY(), negDaughter.dcaZ(), negDaughter.tpcNSigmaPr(), negDaughter.tpcNSigmaPi());
          }
        }
      }
    }
  }

  void init(InitContext&)
  {
    /// init ccdb
    ccdb->setURL(ConfOptions.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    /// init collision selection
    collisionSel.init(ConfCollisionFilter.useEventSel.value);

    /// init track selections
    trackSel.setCheckMinimalSelection(true);
    trackSel.addSelection(ConfTrackBits.tpcClustersMin.value, trackselection::kTPCnClsMin, limits::kLowerLimit, true);
    trackSel.addSelection(ConfTrackBits.tpcCrossedRowsMin.value, trackselection::kTPCcRowsMin, limits::kLowerLimit, true);
    trackSel.addSelection(ConfTrackBits.tpcSharedClustersMax.value, trackselection::kTPCsClsMax, limits::kUpperLimit, true);
    trackSel.addSelection(ConfTrackBits.itsClustersMin.value, trackselection::kITSnClsMin, limits::kLowerLimit, true);
    trackSel.addSelection(ConfTrackBits.itsIbClustersMin.value, trackselection::kITSnClsIbMin, limits::kLowerLimit, true);
    trackSel.addSelection(ConfTrackBits.dcaxyMax.name, ConfTrackFilters.ptMin.value, ConfTrackFilters.ptMax.value, ConfTrackBits.dcaxyMax.value, trackselection::kDCAxyMax, limits::kAbsUpperFunctionLimit, true);
    trackSel.addSelection(ConfTrackBits.dcazMax.name, ConfTrackFilters.ptMin.value, ConfTrackFilters.ptMax.value, ConfTrackBits.dcazMax.value, trackselection::kDCAzMax, limits::kAbsUpperFunctionLimit, true);

    /// init track pid selections
    trackPidSel.setCheckMinimalSelection(false);

    trackPidSel.addSelection(ConfTrackPidBits.itsElectron.value, trackpidselection::kItsElectron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.itsPion.value, trackpidselection::kItsPion, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.itsKaon.value, trackpidselection::kItsKaon, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.itsProton.value, trackpidselection::kItsProton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.itsDeuteron.value, trackpidselection::kItsDeuteron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.itsTriton.value, trackpidselection::kItsTriton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.itsHelium.value, trackpidselection::kItsHelium, limits::kAbsUpperLimit, false);

    trackPidSel.addSelection(ConfTrackPidBits.tpcElectron.value, trackpidselection::kTpcElectron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpcPion.value, trackpidselection::kTpcPion, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpcKaon.value, trackpidselection::kTpcKaon, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpcProton.value, trackpidselection::kTpcProton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpcDeuteron.value, trackpidselection::kTpcDeuteron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpcTriton.value, trackpidselection::kTpcTriton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpcHelium.value, trackpidselection::kTpcHelium, limits::kAbsUpperLimit, false);

    trackPidSel.addSelection(ConfTrackPidBits.tofElectron.value, trackpidselection::kTofElectron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tofPion.value, trackpidselection::kTofPion, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tofKaon.value, trackpidselection::kTofKaon, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tofProton.value, trackpidselection::kTofProton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tofDeuteron.value, trackpidselection::kTofDeuteron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tofTriton.value, trackpidselection::kTofTriton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tofHelium.value, trackpidselection::kTofHelium, limits::kAbsUpperLimit, false);

    trackPidSel.addSelection(ConfTrackPidBits.tpctofElectron.value, trackpidselection::kTpctofElectron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpctofPion.value, trackpidselection::kTpctofPion, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpctofKaon.value, trackpidselection::kTpctofKaon, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpctofProton.value, trackpidselection::kTpctofProton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpctofDeuteron.value, trackpidselection::kTpctofDeuteron, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpctofTriton.value, trackpidselection::kTpctofTriton, limits::kAbsUpperLimit, false);
    trackPidSel.addSelection(ConfTrackPidBits.tpctofHelium.value, trackpidselection::kTpctofHelium, limits::kAbsUpperLimit, false);

    /// init vzero selections
    vzeroSel.setCheckMinimalSelection(true);
    vzeroSel.addSelection(ConfVzeroBits.dcaDaughMax.value, vzeroselection::kDcaDaughMax, limits::kAbsUpperLimit, true);
    vzeroSel.addSelection(ConfVzeroBits.cpaMin.value, vzeroselection::kCpaMin, limits::kLowerLimit, true);
    vzeroSel.addSelection(ConfVzeroBits.transRadMin.value, vzeroselection::kTransRadMin, limits::kLowerLimit, true);
    vzeroSel.addSelection(ConfVzeroBits.transRadMax.value, vzeroselection::kTransRadMax, limits::kUpperLimit, true);
    // vzero positiv daughter selections
    vzeroSel.addSelection(ConfVzeroDaughterBits.dcaMin.value, vzeroselection::kPosDauDcaMin, limits::kLowerLimit, true);
    vzeroSel.addSelection(ConfVzeroDaughterBits.tpcClustersMin.value, vzeroselection::kPosDauTpcClsMin, limits::kLowerLimit, true);
    // vzero negative daughter selections
    vzeroSel.setCheckMinimalSelection(true);
    vzeroSel.addSelection(ConfVzeroDaughterBits.dcaMin.value, vzeroselection::kNegDauDcaMin, limits::kLowerLimit, true);
    vzeroSel.addSelection(ConfVzeroDaughterBits.tpcClustersMin.value, vzeroselection::kNegDauTpcClsMin, limits::kLowerLimit, true);
    // kaon mass limits
    vzeroSel.setKaonMassLimits(ConfVzeroBits.kaonMassRejectionLow.value, ConfVzeroBits.kaonMassRejectionHigh.value);

    /// init vzero daughter pid selections
    vzeroDaugherPidSel.setCheckMinimalSelection(false);
    vzeroDaugherPidSel.addSelection(ConfVzeroDaughterBits.negDaughPionNsigmaMax.value, vzerodaughterpidselection::kNegDaughTpcPion, limits::kAbsUpperLimit, false);
    vzeroDaugherPidSel.addSelection(ConfVzeroDaughterBits.posDaughProtonNsigmaMax.value, vzerodaughterpidselection::kPosDaughTpcProton, limits::kAbsUpperLimit, false);
    vzeroDaugherPidSel.addSelection(ConfVzeroDaughterBits.posDaughPionNsigmaMax.value, vzerodaughterpidselection::kPosDaughTpcPion, limits::kAbsUpperLimit, false);
    vzeroDaugherPidSel.addSelection(ConfVzeroDaughterBits.negDaughProtonNsigmaMax.value, vzerodaughterpidselection::kNegDaughTpcPion, limits::kAbsUpperLimit, false);
  }

  // proccess functions
  // produce tracks for analysis
  void processTracksRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
                           BCsWithTimestamps const&,
                           Filtered<consumeddata::Run3Tracks> const& tracks)
  {
    if (!collisionSel.isSelected<modes::System::kPP_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3Tracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS>(tracksWithItsPid, indexMap);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksRun3pp, "Provide tracks for Run3 analysis", true);

  // produce tracks for analysis (without centrality)
  void proccessTracksRun3ppNoCent(Filtered<consumeddata::Run3PpWithoutCentCollisions>::iterator const& col,
                                  BCsWithTimestamps const&,
                                  Filtered<consumeddata::Run3Tracks> const& tracks)
  {
    if (!collisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3Tracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS>(tracksWithItsPid, indexMap);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, proccessTracksRun3ppNoCent, "Provide tracks for Run3 analysis (use when no centrality calibration is available)", false);

  // produce tracks for QA (without centrality)
  void proccessQaTracksRun3ppNoCent(Filtered<consumeddata::Run3PpWithoutCentCollisions>::iterator const& col,
                                    BCsWithTimestamps const&,
                                    Filtered<consumeddata::Run3TracksFullPid> const& tracks)
  {
    if (!collisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3TracksFullPid, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS_QA>(tracksWithItsPid, indexMap);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, proccessQaTracksRun3ppNoCent, "Provide tracks for Run2 with QA", false);

  // produce tracks and v0s for QA (without centrality)
  void processQaTracksVzerosRun3ppNoCent(Filtered<consumeddata::Run3PpWithoutCentCollisions>::iterator const& col,
                                         BCsWithTimestamps const&,
                                         Filtered<consumeddata::Run3TracksFullPid> const& tracks,
                                         Filtered<consumeddata::Run3PpVzeros> const& v0s)
  {
    if (!collisionSel.isSelected<modes::System::kPP_NoCentCal_Run3>(col)) {
      return;
    }
    initFromCcdb(col.bc_as<BCsWithTimestamps>());

    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = Attach<consumeddata::Run3TracksFullPid, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                   pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);

    fillCollision<modes::System::kPP_NoCentCal_Run3>(col, tracks);
    indexMap.clear();
    fillTracks<modes::Mode::kANALYSIS_QA>(tracksWithItsPid, indexMap);
    fillV0s<modes::Mode::kANALYSIS_QA>(v0s, tracks, indexMap);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processQaTracksVzerosRun3ppNoCent, "Provide Tracks and V0s for Run3 with QA (no centrality calibration)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoUnitedProducer>(cfgc)};
  return workflow;
}
