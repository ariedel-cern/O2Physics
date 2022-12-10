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

/// \file FemtoDreamParticleHisto.h
/// \brief FemtoDreamParticleHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef ANALYSIS_TASKS_PWGCF_O2FEMTODREAM_INCLUDE_O2FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_
#define ANALYSIS_TASKS_PWGCF_O2FEMTODREAM_INCLUDE_O2FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_

#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamParticleHisto
/// \brief Class for histogramming particle properties
/// \tparam particleType Type of the particle (Track/V0/Cascade/...)
/// \tparam suffixType (optional) Takes care of the suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
template <o2::aod::femtodreamparticle::ParticleType particleType, int suffixType = 0>
class FemtoDreamParticleHisto
{
 public:
  /// Destructor
  virtual ~FemtoDreamParticleHisto() = default;

  /// Initialization of the QA histograms
  /// \param registry HistogramRegistry
  template <typename T> 
  void init(HistogramRegistry* registry, T& ptBins, T& DCAxyBins)
  {
    if (registry) {
      mHistogramRegistry = registry;
      
      framework::AxisSpec ptAxis = {ptBins, "#it{p}_{T} (GeV/#it{c})"};
      framework::AxisSpec DCAxyAxis = {DCAxyBins, "DCA_{xy} (cm)"};
      
      /// The folder names are defined by the type of the object and the suffix (if applicable)
      std::string folderName = static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]);
      folderName += static_cast<std::string>(mFolderSuffix[mFolderSuffixType]);

      /// Histograms of the kinematic properties
      mHistogramRegistry->add((folderName + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
      mHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
      mHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});

      /// Particle-type specific histograms
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
        mHistogramRegistry->add((folderName + "/hDCAxy").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
        mHistogramRegistry->add((folderName + "/hCPA").c_str(), "; #it{p}_{T} (GeV/#it{c}); cos#alpha", kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1}});
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }
  
  /// Filling of the histograms
  /// \tparam T Data type of the particle
  /// \param part Particle
  template <typename T>
  void fillQA(T const& part)
  {
    if (mHistogramRegistry) {
      /// Histograms of the kinematic properties
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPt"), part.pt());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hEta"), part.eta());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPhi"), part.phi());

      /// Particle-type specific histograms
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
        mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy"),
                                 part.pt(), part.tempFitVar());
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
        mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hCPA"),
                                 part.pt(), part.tempFitVar());
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }
  
  /// Initialization of the QA histograms for MonteCarlo ()
  /// \param registry HistogramRegistry
  template <typename T> 
  void initMC(HistogramRegistry* registry, T& ptBins, T& DCAxyBins)
  {
    if (registry) {
      mHistogramRegistry = registry;
      
      framework::AxisSpec ptAxis = {ptBins, "#it{p}_{T} (GeV/#it{c})"};
      framework::AxisSpec DCAxyAxis = {DCAxyBins, "DCA_{xy} (cm)"};
      
      /// The folder names are defined by the type of the object and the suffix (if applicable)
      std::string folderName = static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]);
      folderName += static_cast<std::string>(mFolderSuffix[mFolderSuffixType]);
      
      /// Histograms of the kinematic properties
      mHistogramRegistry->add((folderName + "/hPtTruth").c_str(), "; #it{p}_{T} Truth (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
      mHistogramRegistry->add((folderName + "/hEtaTruth").c_str(), "; #eta Truth; Entries", kTH1F, {{200, -1.5, 1.5}});
      mHistogramRegistry->add((folderName + "/hPhiTruth").c_str(), "; #phi Truth; Entries", kTH1F, {{200, 0, 2. * M_PI}});
      
      /// Particle-type specific histograms
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
    
        mHistogramRegistry->add((folderName + "/hDCAxy_Primary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
        mHistogramRegistry->add((folderName + "/hDCAxy_Daughter").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
        mHistogramRegistry->add((folderName + "/hDCAxy_Material").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
        mHistogramRegistry->add((folderName + "/hDCAxy_NotPrimary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
        mHistogramRegistry->add((folderName + "/hDCAxy_Fake").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
        mHistogramRegistry->add((folderName + "/hDCAxy_DaughterLambda").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
        mHistogramRegistry->add((folderName + "/hDCAxy_DaughterSigmaplus").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {ptAxis, DCAxyAxis});
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }
  /// Filling of the histograms for MC
  /// \tparam T Data type of the particle
  /// \param part Particle
  template <typename T>
  void fillQAMC(T const& part)
  {
    if (mHistogramRegistry) {
      
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPtTruth"), part.ptTruth());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hEtaTruth"), part.etaTruth());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPhiTruth"), part.phiTruth());
      
      /// Particle-type specific histograms
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
        switch (part.partOriginMCTruth()){
          case (o2::aod::femtodreamparticleMC::kPrimary):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy_Primary"),
                                part.pt(), part.tempFitVar()); 
            break;
          case (o2::aod::femtodreamparticleMC::kDaughter):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy_Daughter"),
                                part.pt(), part.tempFitVar()); 
            break;
          case (o2::aod::femtodreamparticleMC::kMaterial):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy_Material"),
                                part.pt(), part.tempFitVar()); 
            break;
          case (o2::aod::femtodreamparticleMC::kNotPrimary):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy_NotPrimary"),
                                part.pt(), part.tempFitVar()); 
            break;
          case (o2::aod::femtodreamparticleMC::kFake):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy_Fake"),
                                part.pt(), part.tempFitVar()); 
            break;
          case (o2::aod::femtodreamparticleMC::kDaughterLambda):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy_DaughterLambda"),
                                part.pt(), part.tempFitVar()); 
            break;
          case (o2::aod::femtodreamparticleMC::kDaughterSigmaplus):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy_DaughterSigmaplus"),
                                part.pt(), part.tempFitVar()); 
            break;
          default: 
            LOG(fatal) << "femtodreamparticleMC: not known value for ParticleOriginMCTruth - please check. Quitting!";
        } 
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                   ///< For QA output
  static constexpr o2::aod::femtodreamparticle::ParticleType mParticleType = particleType; ///< Type of the particle under analysis
  static constexpr int mFolderSuffixType = suffixType;                                     ///< Counter for the folder suffix specified below
  static constexpr std::string_view mFolderSuffix[5] = {"", "_one", "_two", "_oneMC", "_twoMC"};               ///< Suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
};
} // namespace o2::analysis::femtoDream

#endif /* ANALYSIS_TASKS_PWGCF_O2FEMTODREAM_INCLUDE_O2FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_ */