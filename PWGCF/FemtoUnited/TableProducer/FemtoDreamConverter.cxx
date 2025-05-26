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

/// \file FemtoDreamConverter.cxx
/// \brief converter task for femtoDream to femtoUnited
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoVzerosDerived.h"

#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FemtoDreamConverter {

  // femto collisions
  Produces<FUCols> outputCollisions;

  // femto tracks
  Produces<FUTracks> outputTracks;
  Produces<FUTrackMasks> outputTrackMasks;
  // Produces<FTrackDCAs> outputTrackDCAs;
  // Produces<FTrackExtras> outputTrackExtras;
  //
  // // femto vzeros
  Produces<FUVzeros> outputVzeros;
  // Produces<FVzeroExtras> outputVzeroExtras;
  //
  Produces<FUVzeroDaus> outputVzeroDaus;

  // void init(o2::framework::InitContext&){};

  void process(FDCollision const& col, FDParticles const& parts)
  {
    outputCollisions(col.posZ(), col.multNtr(), col.multV0M(), col.sphericity(), col.magField());
    for (auto const& part : parts) {
      // fill tracks
      if (part.partType() == femtodreamparticle::ParticleType::kTrack) {
        // use last index for the collision index
        // the index stored in old format might be different, but it is only about having a
        // unique index for every collision
        bool hasPositiveCharge = false;
        if ((part.cut() & 0b11) == 0b10) {
          hasPositiveCharge = true;
        }
        outputTracks(outputCollisions.lastIndex(), part.pt(), part.eta(), part.phi(), hasPositiveCharge);
        // in old format tpc and tpctof information is stored in the same bitmask
        // so we only fill the TPC mask in the FemtoUnited framework
        // for analysis the bitmask has to be computed with the old tools
        outputTrackMasks(part.cut(), part.pidcut());

        // fill Vzeros
        // when Vzeros are filled the daughters will also be filled, so they are skipped in the outer loop
      } else if (part.partType() == femtodreamparticle::ParticleType::kV0) {
        // vzero
        outputVzeros(outputCollisions.lastIndex(), part.pt(), part.eta(), part.phi(), part.mLambda(), part.mAntiLambda());
        // vzero daughters
        const auto& posChild = parts.iteratorAt(part.index() - 2);
        const auto& negChild = parts.iteratorAt(part.index() - 1);
        outputVzeroDaus(
          posChild.childrenIds()[0], posChild.pt(), posChild.eta(), posChild.phi(),
          negChild.childrenIds()[1], negChild.pt(), negChild.eta(), negChild.phi());
      }
    }
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoDreamConverter>(cfgc)};
  return workflow;
}
