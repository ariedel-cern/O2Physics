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

#ifndef PWGCF_FEMTOUNITED_CORE_MODES_H_
#define PWGCF_FEMTOUNITED_CORE_MODES_H_

#include <cstdint>

namespace femtomodes
{

enum class Mode : uint32_t {
  kANALYSIS = 0x1,
  kQA = 0x2,
  kMC = 0x4,
  kANALYSIS_QA = kANALYSIS | kQA,
  kANALYSIS_MC = kANALYSIS | kMC,
  kANALYSIS_QA_MC = kANALYSIS | kQA | kMC,
};

// Function to check if a mode is activated
constexpr bool isModeSet(Mode mode, Mode flag)
{
  return static_cast<uint32_t>(mode) & static_cast<uint32_t>(flag);
}

enum class System : uint32_t {
  kPP = 0x1,
  kPbPb = 0x2
};

// Function to check if a mode is activated
constexpr bool isSystemSet(System mode, System flag)
{
  return static_cast<uint32_t>(mode) & static_cast<uint32_t>(flag);
}

} // namespace femtomodes

#endif // PWGCF_FEMTOUNITED_CORE_MODES_H_
