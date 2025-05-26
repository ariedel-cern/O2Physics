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

#include <cmath>
#include <bitset>
#include <vector>
#include <algorithm>
#include <string>

#include "TF1.h"

#include "fairlogger/Logger.h"
#include "CommonConstants/MathConstants.h"

namespace o2::analysis::femtounited
{

/// Limit type for selections
namespace limits
{
enum LimitType { kUpperLimit,            ///< simple upper limit for the value, e.g. p_T < 1 GeV/c
                 kAbsUpperLimit,         ///< upper limit of the absolute value, e.g. |eta| < 0.8
                 kLowerLimit,            ///< simple lower limit for the value, e.g. p_T > 0.2 GeV/c
                 kAbsLowerLimit,         ///< lower limit of the absolute value, e.g. |DCA_xyz| > 0.05 cm
                 kEqual,                 ///< values need to be equal, e.g. sign = 1
                 kUpperFunctionLimit,    ///< simple upper limit of a function value, e.g. DCA_xy > f(pt)
                 kAbsUpperFunctionLimit, ///< upper limit of an absolute value given by a function, e.g. |DCA_xy| > f(pt)
                 kLowerFunctionLimit,    ///< simple upper limit of a function value, e.g. DCA_xy > f(pt)
                 kAbsLowerFunctionLimit  ///< upper limit of an absolute value given by a function, e.g. |DCA_xy| > f(pt)
};
} // namespace limits

/// Simple class for storing selections of a single observable
/// \tparam T Data type used for the selection values (float/int/...)
/// \tparam BitmaskType Compute number of selections from BitmaskType (should be some unsigned integer)
template <typename T, typename BitmaskType>
class SelectionContainer
{
 public:
  /// Default constructor
  SelectionContainer() {}

  /// Constructor
  /// \param values Vector of values for the selection
  /// \param limitType Type of limit of the selection
  /// \param SkipLastBit Boolean whether to skip the last bit
  SelectionContainer(std::vector<T>& values, limits::LimitType limitType, bool SkipLastBit)
    : mValues(values),
      mLimitType(limitType),
      mSkipLastBit(SkipLastBit)
  {
    if (mValues.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections for single a observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    // if limitTypeBoolean whether to skip the last bits kEqual we can never skip the last bit
    if (limitType == limits::kEqual) {
      mSkipLastBit = false;
    }
    // values for selection are not necessarily ordered correctly
    sortSelections();
  }

  /// Constructor
  /// \param baseName base name for TF1 object
  /// \param lowerLimit upper limit of the TF1 object
  /// \param upperLimit lower limit of the TF1 object
  /// \param function vector of strings for initializing of the TF1 object
  /// \param limitType Type of limit of the selection
  /// \param SkipLastBit Boolean whether to skip the last bit
  SelectionContainer(std::string baseName, T lowerLimit, T upperLimit, std::vector<std::string>& functions, limits::LimitType limitType, bool SkipLastBit)
    : mLimitType(limitType), mSkipLastBit(SkipLastBit)
  {
    if (mValues.size() > sizeof(BitmaskType) * CHAR_BIT) {
      LOG(fatal) << "Too many selections for single a observable. Limit is " << sizeof(BitmaskType) * CHAR_BIT;
    }
    for (std::size_t i = 0; i < functions.size(); i++) {
      mFunctions.emplace_back((baseName + std::to_string(i)).c_str(), functions.at(i).c_str(), lowerLimit, upperLimit);
    }
    // functions for selection are not necessarily ordered correctly
    // use value at midpoint to order them
    // here we rely on the user that the functions can be ordered like this over the whole interval
    T midPoint = (lowerLimit + upperLimit) / 2.;
    sortFunctions(midPoint);
    // initialize the values also to the midpoint
    this->updateLimits(midPoint);
  }

  /// Destructor
  virtual ~SelectionContainer() = default;

  /// Sort selections accroding to limit type
  void sortSelections()
  {
    switch (mLimitType) {
      case (limits::kUpperLimit):
      case (limits::kAbsUpperLimit):
        std::sort(mValues.begin(), mValues.end(), [](T a, T b) { return a >= b; });
        break;
      case (limits::kLowerLimit):
      case (limits::kAbsLowerLimit):
      case (limits::kEqual):
        std::sort(mValues.begin(), mValues.end(), [](T a, T b) { return a <= b; });
        break;
      default:
        break;
    }
  }

  // sort limit functions
  void sortFunctions(T value)
  {
    switch (mLimitType) {
      case (limits::kUpperFunctionLimit):
      case (limits::kAbsUpperFunctionLimit):
        std::sort(mFunctions.begin(), mFunctions.end(), [value](TF1 a, TF1 b) { return a.Eval(value) >= b.Eval(value); });
        break;
      case (limits::kLowerFunctionLimit):
      case (limits::kAbsLowerFunctionLimit):
        std::sort(mFunctions.begin(), mFunctions.end(), [value](TF1 a, TF1 b) { return a.Eval(value) <= b.Eval(value); });
        break;
      default:
        break;
    }
  }

  // update the selection limits depending on the passed function
  void updateLimits(T value)
  {
    // functions are ordered so just add the values in the same order
    for (std::size_t i = 0; i < mValues.size(); i++) {
      mValues.at(i) = mFunctions.at(i).Eval(value);
    }
  }

  /// Check which selections are fulfilled
  /// \param observable Value of the variable to be checked
  void setBitmask(T variable)
  {
    // better safe than sorry and reset the bitmask before you evaluate a new observable
    mBitmask.reset();
    // the values are ordered, for most loost to  most tight, as soon as one comparison is not true, we can break out of the loop
    bool breakLoop = false;
    // iterate over all limits and set the corresponding bit if we pass the selection, otherwise break out as soon as we can
    for (size_t i = 0; i < mValues.size(); i++) {
      switch (mLimitType) {
        case (limits::kUpperLimit):
        case (limits::kUpperFunctionLimit):
          if (variable <= mValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kAbsUpperLimit):
        case (limits::kAbsUpperFunctionLimit):
          if (std::abs(variable) <= mValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kLowerLimit):
        case (limits::kLowerFunctionLimit):
          if (variable >= mValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kAbsLowerLimit):
        case (limits::kAbsLowerFunctionLimit):
          if (std::abs(variable) >= mValues.at(i)) {
            mBitmask.set(i);
          } else {
            breakLoop = true;
          }
          break;
        case (limits::kEqual):
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
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> getBitmask() const
  {
    // if we do not skip the last bit, return full bitmask
    if (mSkipLastBit == false) {
      return mBitmask;
    }
    // if we do, we also need to check for kEqual first
    if (limits::kEqual == mLimitType) {
      // if the limit is equal, return all the bits
      return mBitmask;
    } else {
      // for the other selections we can remove the first bit since it is the minimal selection and therefore always true
      return mBitmask >> 1;
    }
  }

  /// Check whether the minimal selection is fulfilled or not
  /// \return Whether the selection is fulfilled or not
  bool getMinimalSelection() const { return mBitmask.any(); }

  /// Return loosest selection value
  ///  The values are ordered, the loosest value will be the first one
  /// \return loosest selection
  T getLoosestSelection() const { return mValues.at(0); }

  /// Check whether selections are empty
  bool empty() const { return mValues.empty(); }

  /// Get shift for final bitmask
  int getShift() const
  {
    if (mValues.empty()) {
      return 0;
    }
    if (mSkipLastBit) {
      return static_cast<int>(mValues.size() - 1);
    } else {
      return static_cast<int>(mValues.size());
    }
  }

 private:
  std::vector<T> mValues{};                             ///< Values used for the selection
  std::vector<TF1> mFunctions{};                        ///< Values used for the selection
  limits::LimitType mLimitType;                         ///< Limit type of selection
  std::bitset<sizeof(BitmaskType) * CHAR_BIT> mBitmask; ///< bitmask for a given observable
  bool mSkipLastBit = false;
};

} // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_SELECTIONCONTAINER_H_
