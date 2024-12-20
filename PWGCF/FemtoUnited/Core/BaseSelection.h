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

/// \file FemtoDreamObjectSelection.h
/// \brief FemtoDreamObjectSelection - Parent class of all selections
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_

#include <algorithm>
#include <vector>

#include "fairlogger/Logger.h"

#include "PWGCF/FemtoUnited/Core/SelectionContainer.h"

namespace o2::analysis::femtounited
{
/// \class BaseSelection
/// \brief Base class to contain all cuts and assemble bitmask
/// \tparam valueType Data type used for the selection (float/int/...)
/// \tparam BitmaskSize Size of the bitmask
template <typename T, typename BitmaskType, size_t NObservables>
class BaseSelection
{
 public:
  /// Constructor
  BaseSelection() {}

  /// Destructor
  virtual ~BaseSelection() = default;

  /// Pass the Configurable of selection values in the analysis task to the selection class
  /// \param configSelection Vector from configurable containing the values employed for the selection
  /// \param observableType Observable to be employed for the selection
  /// \param limitType Type of the selection limit
  void addSelection(std::vector<T>& configSelections, int Observable, limits::LimitType limitType, bool SkipLastBit)
  {
    if (static_cast<size_t>(Observable) >= NObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NObservables;
    }
    mNSelections += configSelections.size();
    if (mNSelections >= 8 * sizeof(BitmaskType)) {
      LOG(fatal) << "Too many selections. At most " << 8 * sizeof(BitmaskType) << " are supported";
    }
    mSelections.at(Observable) = SelectionContainer<T>(configSelections, limitType, SkipLastBit);
  }

  void resetMinimalSelection() { mMinimalSelected = true; }

  void setCheckMinimalSelection(bool checkMinimalSelection) { mCheckMinimalSelection = checkMinimalSelection; }

  /// set bitmask for a given observable
  /// \param observable Observable to be checked
  /// \param value Value of the observable
  void setBitmaskForObservable(int observable, T value)
  {
    // if any object did not pass minimal selections, there is no point in setting bitmask for other observables
    // minimal selection for each observable is computed after adding it
    // can be deactivate by setting mCheckMinimalSelection to false
    if (mCheckMinimalSelection == true && mMinimalSelected == false) {
      return;
    }

    // set bitmask for given observable
    mSelections.at(observable).setBitmask(value);
    // check if minimal selction for this observable holds
    if (mSelections.at(observable).getMinimalSelection() == false) {
      mMinimalSelected = false;
    }
  }

  /// check if minimal Selections are passed
  bool getMinimalSelection() { return mMinimalSelected; }

  /// check if any Selections are passed
  bool getAnySelection()
  {
    return std::any_of(mSelections.begin(), mSelections.end(), [](auto& selection) {
      return selection.getMinimalSelection();
    });
  }

  /// assemble final bitmask
  void assembleBismask()
  {
    // if minimal selections are not passed, just set bitmask to 0
    if (mCheckMinimalSelection == true && mMinimalSelected == false) {
      mBitmask.reset();
    }

    // to assemble bitmask, convert all bitmask into integers
    uint64_t result = 0u;
    int shift = 0;
    uint64_t value = 0u;
    for (auto& selection : mSelections) {
      // if there are no selections for a certain observable, skip
      if (selection.empty()) {
        continue;
      }
      // Convert the bitset to an integer
      value = selection.getBitmask().to_ullong();
      shift = selection.getShift();
      // Shift the result to make space and add the new value
      result = (result << shift) | value;
    }
    // Convert the accumulated integer result to a bitset
    mBitmask = std::bitset<8 * sizeof(BitmaskType)>(result);
  }

  BitmaskType getBitmask() { return static_cast<BitmaskType>(mBitmask.to_ullong()); }

 protected:
  std::array<SelectionContainer<T>, NObservables> mSelections; ///< Array containing all selections
  std::bitset<8 * sizeof(BitmaskType)> mBitmask;               ///< final bitmaks
  size_t mNSelections = 0;                                     ///< Number of selections
  bool mMinimalSelected = true;                                ///< whether minimal selections are fullfilled
  bool mCheckMinimalSelection = false;                         ///< whether to check for minimal selections bitmaks
};

} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_
