// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUtils.h
/// \brief Collision selection
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_UTILS_FEMTOUTILS_H_
#define PWGCF_FEMTOUNITED_UTILS_FEMTOUTILS_H_

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>
#include "TPDGCode.h"
#include "CommonConstants/PhysicsConstants.h"

namespace o2::analysis::femtounited
{
namespace utils
{

inline float geometricMean(float a, float b)
{
  return std::sqrt(a * a + b * b);
}

inline float dca(float dcaXY, float dcaZ)
{
  return std::sqrt(dcaXY * dcaXY + dcaZ * dcaZ);
}

template <typename T>
T getDaughterIndex(T daugherIndex, std::vector<std::pair<T, T>> map)
{
  for (const auto& pair : map) {
    if (pair.first == daugherIndex) {
      return pair.second; // Return the second element of the found pair
    }
  }
  return -1;
};

template <typename T>
float itsSignal(T const& track)
{
  uint32_t clsizeflag = track.itsClusterSizes();
  auto clSizeLayer0 = (clsizeflag >> (0 * 4)) & 0xf;
  auto clSizeLayer1 = (clsizeflag >> (1 * 4)) & 0xf;
  auto clSizeLayer2 = (clsizeflag >> (2 * 4)) & 0xf;
  auto clSizeLayer3 = (clsizeflag >> (3 * 4)) & 0xf;
  auto clSizeLayer4 = (clsizeflag >> (4 * 4)) & 0xf;
  auto clSizeLayer5 = (clsizeflag >> (5 * 4)) & 0xf;
  auto clSizeLayer6 = (clsizeflag >> (6 * 4)) & 0xf;
  int numLayers = 7;
  int sumClusterSizes = clSizeLayer1 + clSizeLayer2 + clSizeLayer3 + clSizeLayer4 + clSizeLayer5 + clSizeLayer6 + clSizeLayer0;
  float cosLamnda = 1. / std::cosh(track.eta());
  return (static_cast<float>(sumClusterSizes) / numLayers) * cosLamnda;
};

template <typename T>
float sphericity(T const& tracks)
{
  if (tracks.size() <= 2) {
    return 2.;
  }

  // Initialize the transverse momentum tensor components
  float sxx = 0.;
  float syy = 0.;
  float sxy = 0.;
  float sumPt = 0.;

  // Loop over the tracks to compute the tensor components
  for (const auto& track : tracks) {
    sxx += (track.px() * track.px()) / track.pt();
    syy += (track.py() * track.py()) / track.pt();
    sxy += (track.px() * track.py()) / track.pt();
    sumPt += track.pt();
  }
  sxx /= sumPt;
  syy /= sumPt;
  sxy /= sumPt;

  // Compute the eigenvalues (real values)
  float lambda1 = ((sxx + syy) + std::sqrt((sxx + syy) * (sxx + syy) - 4 * (sxx * syy - sxy * sxy))) / 2;
  float lambda2 = ((sxx + syy) - std::sqrt((sxx + syy) * (sxx + syy) - 4 * (sxx * syy - sxy * sxy))) / 2;

  if (lambda1 <= 0 || lambda2 <= 0) {
    return 2.;
  }

  // Compute sphericity
  return 2. * lambda2 / (lambda1 + lambda2);
}

inline float getMass(int pdgCode)
{
  // use this function instead of TDatabasePDG to return masses defined in the PhysicsConstants.h header
  // this approach saves a lot of memory and important partilces like deuteron are missing in TDatabasePDG anyway
  float mass = 0;
  // add new particles if necessary here
  switch (std::abs(pdgCode)) {
    case kPiPlus:
      mass = o2::constants::physics::MassPiPlus;
      break;
    case kKPlus:
      mass = o2::constants::physics::MassKPlus;
      break;
    case kProton:
      mass = o2::constants::physics::MassProton;
      break;
    case kLambda0:
      mass = o2::constants::physics::MassLambda;
      break;
    case o2::constants::physics::Pdg::kPhi:
      mass = o2::constants::physics::MassPhi;
      break;
    case o2::constants::physics::Pdg::kLambdaCPlus:
      mass = o2::constants::physics::MassLambdaCPlus;
      break;
    case o2::constants::physics::Pdg::kDeuteron:
      mass = o2::constants::physics::MassDeuteron;
      break;
    case o2::constants::physics::Pdg::kTriton:
      mass = o2::constants::physics::MassTriton;
      break;
    case o2::constants::physics::Pdg::kHelium3:
      mass = o2::constants::physics::MassHelium3;
      break;
    default:
      LOG(fatal) << "PDG code is not suppored";
  }
  return mass;
}

inline float dphistar(float magfield, float radius, float charge, float pt, float phi)
{
  float arg = 0.3 * charge * magfield * radius * 0.01 / (2. * pt);
  // for very low pT particles, this value goes outside of range -1 to 1 at at large tpc radius; asin fails
  if (std::abs(arg) < 1.f) {
    return phi - std::asin(arg);
  } else {
    return 99.f;
  }
}

}; // namespace utils
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_UTILS_FEMTOUTILS_H_
