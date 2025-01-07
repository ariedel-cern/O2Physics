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

/// \file ProducerHelpers.h
/// \brief Collision selection
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_UTLITIES_FEMTOUTILS_H_
#define PWGCF_FEMTOUNITED_UTLITIES_FEMTOUTILS_H_

#include <cmath>
#include "CommonConstants/MathConstants.h"

#include "CCDB/BasicCCDBManager.h"

#include "PWGCF/FemtoUnited/Core/Modes.h"

namespace o2::analysis::femtounited
{
namespace utils
{

inline float geometricMean(float a, float b)
{
  return std::sqrt(a * a + b * b);
}

inline float Dca(float dcaXY, float dcaZ)
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
  float cos_lambda = 1. / std::cosh(track.eta());
  return (static_cast<float>(sumClusterSizes) / numLayers) * cos_lambda;
};

template <typename T>
float sphericity(T const& tracks)
{
  if (tracks.size() <= 2) {
    return 2.;
  }

  // Initialize the transverse momentum tensor components
  float Sxx = 0.;
  float Syy = 0.;
  float Sxy = 0.;
  float SumPt = 0.;

  // Loop over the tracks to compute the tensor components
  for (const auto& track : tracks) {
    Sxx += (track.px() * track.px()) / track.pt();
    Syy += (track.py() * track.py()) / track.pt();
    Sxy += (track.px() * track.py()) / track.pt();
    SumPt += track.pt();
  }
  Sxx /= SumPt;
  Syy /= SumPt;
  Sxy /= SumPt;

  // Compute the eigenvalues (real values)
  float lambda1 = ((Sxx + Syy) + std::sqrt((Sxx + Syy) * (Sxx + Syy) - 4 * (Sxx * Syy - Sxy * Sxy))) / 2;
  float lambda2 = ((Sxx + Syy) - std::sqrt((Sxx + Syy) * (Sxx + Syy) - 4 * (Sxx * Syy - Sxy * Sxy))) / 2;

  if (lambda1 <= 0 || lambda2 <= 0) {
    return 2.;
  };

  // Compute sphericity
  return 2. * lambda2 / (lambda1 + lambda2);
}

}; // namespace utils
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_COLSELECTION_H_
