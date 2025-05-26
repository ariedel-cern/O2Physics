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

/// \file BaseSelection.h
/// \brief Definition of the BaseSelection class
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_

#include <algorithm>
#include <string>
#include <vector>

#include "fairlogger/Logger.h"

#include "PWGCF/FemtoUnited/Core/SelectionContainer.h"

namespace o2::analysis::femtounited
{
/// \class BaseSelection
/// \brief Base class to contain all cuts and assemble bitmask
/// \tparam T Data type used for the selections (float/int/...)
/// \tparam BitmaskType Compute size of the bitmask from BitmaskType
/// \tparam NObservables Number of observables
template <typename T, typename BitmaskType, size_t NObservables>
class BaseSelection
{
 public:
  /// Constructor
  BaseSelection() {}

  /// Destructor
  virtual ~BaseSelection() = default;

  /// Pass the Configurable of selection values in the analysis task to the selection class
  /// \param configSelections Vector of configurables containing the values employed for the selection
  /// \param Observable Observable to be employed for the selection
  /// \param limitType Type of the selection limit
  void addSelection(std::vector<T>& configSelections, int Observable, limits::LimitType limitType, bool SkipLastBit)
  {
    if (static_cast<size_t>(Observable) >= NObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NObservables;
    }
    mNSelections += configSelections.size();
    if (mNSelections >= sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
    mSelections.at(Observable) = SelectionContainer<T, BitmaskType>(configSelections, limitType, SkipLastBit);
  }

  /// Pass the Configurable of selection values in the analysis task to the selection class
  /// \param configSelection Vector from configurable containing the values employed for the selection
  /// \param observableType Observable to be employed for the selection
  /// \param limitType Type of the selection limit
  void addSelection(std::string baseName, T lowerLimit, T upperLimit, std::vector<std::string>& configSelections, int Observable, limits::LimitType limitType, bool SkipLastBit)
  {
    if (static_cast<size_t>(Observable) >= NObservables) {
      LOG(fatal) << "Observable is not valid. Observable (index) has to be smaller than " << NObservables;
    }
    mNSelections += configSelections.size();
    if (mNSelections >= sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections. At most " << sizeof(BitmaskType) * CHAR_BIT << " are supported";
    }
    mSelections.at(Observable) = SelectionContainer<T, BitmaskType>(baseName, lowerLimit, upperLimit, configSelections, limitType, SkipLastBit);
  }

  void updateLimits(int observable, T value) { mSelections.at(observable).updateLimits(value); }

  void reset()
  {
    mBitmask.reset();
    mMinimalSelected = true;
  }

  void setCheckMinimalSelection(bool checkMinimalSelection) { mCheckMinimalSelection = checkMinimalSelection; }

  /// set bitmask for a given observable
  /// \param observable Observable to be checked
  /// \param value Value of the observable
  void setBitmaskForObservable(int observable, T value)
  {
    // if there are no values set, bail out
    if (mSelections.at(observable).empty()) {
      return;
    }
    // if any previous observable did not pass minimal selections, there is no point in setting bitmask for other observables
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

  bool getAnySelection(int Observable)
  {
    return mSelections.at(Observable).getMinimalSelection();
  }

  /// assemble final bitmask
  void assembleBismask()
  {
    // if minimal selections are not passed, just set bitmask to 0
    if (mCheckMinimalSelection == true && mMinimalSelected == false) {
      mBitmask.reset();
      return;
    }

    // to assemble bitmask, convert all bitmask into integers
    // shift the current one and add the new bits
    // uint64_t result = 0u;
    // int shift = 0;
    // uint64_t value = 0u;
    for (auto const& selection : mSelections) {
      // if there are no selections for a certain observable, skip
      if (selection.empty()) {
        continue;
      }
      // Convert the bitset to an integer
      // value = selection.getBitmask().to_ullong();
      // shift = selection.getShift();
      // Shift the result to make space and add the new value
      mBitmask = (mBitmask << selection.getShift()) | selection.getBitmask();
    }
    // Convert the accumulated integer result to a bitset
    // mBitmask = std::bitset<CHAR_BIT * sizeof(BitmaskType)>(result);
  }

  BitmaskType getBitmask() { return static_cast<BitmaskType>(mBitmask.to_ullong()); }

 protected:
  std::array<SelectionContainer<T, BitmaskType>, NObservables> mSelections; ///< Array containing all selections
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mBitmask;                     ///< final bitmaks
  size_t mNSelections = 0;                                                  ///< Number of selections
  bool mMinimalSelected = true;                                             ///< whether minimal selections are fullfilled
  bool mCheckMinimalSelection = false;                                      ///< whether to check for minimal selections bitmaks
};

} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_BASESELECTION_H_
