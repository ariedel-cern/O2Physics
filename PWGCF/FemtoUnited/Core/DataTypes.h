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

#ifndef PWGCF_FEMTOUNITED_CORE_DATATYPES_H_
#define PWGCF_FEMTOUNITED_CORE_DATATYPES_H_

#include <cstdint>

namespace o2::aod
{
namespace femtodatatypes
{
// Note: Length of the bitmask is the limit of how many selections can be configured

// Bitmaks for tracks
using TrackMaskType = uint32_t;
using TrackPidMaskType = uint32_t;
// using TrackTPCMaskType = uint16_t;
// using TrackTOFMaskType = uint16_t;
// using TrackTPCTOFMaskType = uint16_t;

// Bitmaks for vzeros and daughters
using VzeroMaskType = uint32_t;
// using VzeroDauTrackMaskType = uint8_t;
// using VzeroDauTPCMaskType = uint8_t;

// Bitmaks for cascades, vzero daughter and bachelor
using CascadeMaskType = uint32_t;
// using CascadeBachelorMaskType = uint16_t;
// using CascadeBachelorTPCMaskType = uint16_t;
// using VzeroDauTrackMaskType = uint8_t;
// using VzeroDauTPCMaskType = uint8_t;

} // namespace femtodatatypes

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_CORE_DATATYPES_H_
