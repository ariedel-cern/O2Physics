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

/// \file Histmanager.h
/// \brief Femtounited collision histograms

#ifndef PWGCF_FEMTOUNITED_CORE_HISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_HISTMANAGER_H_

#include "Framework/HistogramRegistry.h"

namespace Histmanager
{

template <typename Hist>
struct HistInfo {
  Hist hist;
  o2::framework::HistType histtype;
  std::string_view histname;
  std::string_view histdesc;
};

template <typename EnumType, typename ArrayType>
constexpr o2::framework::HistType GetHistType(EnumType variable, const ArrayType& array)
{
  for (const auto& entry : array) {
    if (entry.hist == variable) {
      return entry.histtype;
    }
  }
  return o2::framework::kUndefinedHist;
}

template <typename EnumType, typename ArrayType>
constexpr std::string_view GetHistName(EnumType variable, const ArrayType& array)
{
  for (const auto& entry : array) {
    if (entry.hist == variable) {
      return entry.histname;
    }
  }
  return ""; // Return an empty string or a default name if not found
}

template <typename EnumType, typename ArrayType>
constexpr std::string GetHistNamev2(EnumType variable, const ArrayType& array)
{
  for (const auto& entry : array) {
    if (entry.hist == variable) {
      return std::string(entry.histname);
    }
  }
  return std::string(""); // Return an empty string or a default name if not found
}

template <typename EnumType, typename ArrayType>
constexpr const char* GetHistDesc(EnumType variable, const ArrayType& array)
{
  for (const auto& entry : array) {
    if (entry.hist == variable) {
      return entry.histdesc.data();
    }
  }
  return ""; // Return an empty string or a default description if not found
}

} // namespace Histmanager

#endif
