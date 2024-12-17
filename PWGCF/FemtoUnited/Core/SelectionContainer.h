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

/// \file SelectionContainer.h
/// \brief SelectionContainer - small container holding selections of an observable
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_SELECTIONCONTAINER_H_
#define PWGCF_FEMTOUNITED_CORE_SELECTIONCONTAINER_H_

#include "fairlogger/Logger.h"
#include <cmath>
#include <bitset>
#include <vector>
#include <algorithm>
#include "CommonConstants/MathConstants.h"

namespace o2::analysis::femtounited
{

namespace Limits
{
/// Limit type for selections
enum LimitType { kUpperLimit,    ///< simple upper limit for the value, e.g. p_T < 1 GeV/c
                 kAbsUpperLimit, ///< upper limit of the absolute value, e.g. |eta| < 0.8
                 kLowerLimit,    ///< simple lower limit for the value, e.g. p_T > 0.2 GeV/c
                 kAbsLowerLimit, ///< lower limit of the absolute value, e.g. |DCA_xyz| > 0.05 cm
                 // kUpperFunctionLimit,    ///< simple upper limit of a function value, e.g. DCA_xy > f(pt)
                 // kAbsUpperFunctionLimit, ///< upper limit of an absolute value given by a function, e.g. |DCA_xy| > f(pt)
                 kEqual ///< values need to be equal, e.g. sign = 1
};

} // namespace Limits

// bitsets need number of bits at compile time. Set reasonable limit here
// This limits the number of selections for ONE observable
constexpr size_t BitmaskMaxSize = 16;

/// Simple class for storing selections of a single observable
/// \tparam SelDataType Data type used for the selection values (float/int/...)
template <typename T>
class SelectionContainer
{
 public:
  /// Default constructor
  SelectionContainer() {};

  /// Constructor
  /// \param values Values for the selection
  /// \param limitType Type of limit of the selection
  SelectionContainer(std::vector<T> values, Limits::LimitType limitType, bool SkipLastBit)
    : mValues(values),
      mLimitType(limitType),
      mSkipLastBit(SkipLastBit)
  {
    if (mValues.size() > BitmaskMaxSize) {
      LOG(fatal) << "Too many selections for single a observable. Current limit is " << BitmaskMaxSize;
    }
    if (limitType == Limits::kEqual) {
      mSkipLastBit = false;
    }
    // values for selection are not necessarily ordered correctly
    sortSelections();
  }

  /// Destructor
  virtual ~SelectionContainer() = default;

  /// Sort selections accroding to limit type
  void sortSelections()
  {
    switch (mLimitType) {
      case (Limits::LimitType::kUpperLimit):
      case (Limits::LimitType::kAbsUpperLimit):
        std::sort(mValues.begin(), mValues.end(), [](T a, T b) { return a >= b; });
        break;
      case (Limits::LimitType::kLowerLimit):
      case (Limits::LimitType::kAbsLowerLimit):
      case (Limits::LimitType::kEqual):
        std::sort(mValues.begin(), mValues.end(), [](T a, T b) { return a <= b; });
        break;
    }
  }

  /// Check which selections are fulfilled
  /// \param observable Value of the variable to be checked
  void setBitmask(T variable)
  {
    // better safe than sorry and reset the bitmask before you evaluate a new observable
    mBitmask.reset();
    // the values are order, as soon as one comparison is not true, we can break out of the loop
    bool breakLoop = false;
    // iterate over all limits and set the corresponding bit if we pass the selection, otherwise break out as soon as we can
    for (size_t i = 0; i < mValues.size(); i++) {
      switch (mLimitType) {
        case (Limits::LimitType::kUpperLimit):
          if (variable <= mValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (Limits::LimitType::kAbsUpperLimit):
          if (std::abs(variable) <= mValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (Limits::LimitType::kLowerLimit):
          if (variable >= mValues.at(i)) {
            mBitmask.set(i);
          } else
            breakLoop = true;
          break;
        case (Limits::LimitType::kAbsLowerLimit):
          if (std::abs(variable) >= mValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (Limits::LimitType::kEqual):
          // special case for kEqual since here we cannot really establish an order so we need to check all cases explicitly
          if (std::abs(variable - mValues.at(i)) < constants::math::Epsilon) {
            mBitmask.set(i);
          }
          break;
        default:
          breakLoop = true;
      }
      if (breakLoop) {
        break;
      }
    }
  }

  /// Return the bitmask of the selections
  /// \return bitmask
  std::bitset<BitmaskMaxSize> getBitmask()
  {
    // if we do not skip the last bit, return full bitmask
    if (mSkipLastBit == false) {
      return mBitmask;
    }

    // if we do, we also need to check for kEqual first
    if (Limits::kEqual == mLimitType) {
      // if the limit is equal, return all the bits
      return mBitmask;
    } else {
      // for the other selections we can remove the first bit since it is the minimal selection and therefore always true
      return mBitmask >> 1;
    }
  }

  /// Check whether the minimal selection is fulfilled or not
  /// \return Whether the selection is fulfilled or not
  bool getMinimalSelection() { return mBitmask.any(); }

  /// Check whether selections are empty
  bool empty() { return mValues.empty(); }

  /// Get shift for final bitmask
  int getShift()
  {
    if (mSkipLastBit) {
      return static_cast<int>(mValues.size() - 1);
    } else {
      return static_cast<int>(mValues.size());
    }
  }

 private:
  std::vector<T> mValues{};             ///< Values used for the selection
  Limits::LimitType mLimitType;         ///< Limit type of selection
  std::bitset<BitmaskMaxSize> mBitmask; ///< bitmask for a given observable
  bool mSkipLastBit = false;
};

} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_SELECTIONCONTAINER_H_
