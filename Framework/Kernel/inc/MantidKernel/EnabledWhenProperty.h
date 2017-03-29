#ifndef MANTID_KERNEL_ENABLEDWHENPROPERTY_H_
#define MANTID_KERNEL_ENABLEDWHENPROPERTY_H_

#include "MantidKernel/System.h"
#include "MantidKernel/IPropertyManager.h"
#include "MantidKernel/IPropertySettings.h"
#include <memory>

namespace Mantid {
namespace Kernel {

/** IPropertySettings for a property that sets it to enabled (in the GUI)
when the value of another property is:
- its default (or not)
- equal to a string (or not)

Usage:

- In an algorithm's init() method, after a call to create a property:

declareProperty("PropA", 123);

- Add a call like this:

setPropertySettings("PropA",
make_unique<EnabledWhenProperty>("OtherProperty",
IS_EQUAL_TO, "2000");

- This will make the property "PropA" show as enabled when
"OtherProperty"'s value is equal to "2000". Similarly, you can use:

setPropertySettings("PropA",
make_unique<VisibleWhenProperty>("OtherProperty",
IS_NOT_DEFAULT);

- This will make the property "PropA" show as visible when "OtherProperty"
is NOT the default value for it.


Copyright &copy; 2011-2017 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
National Laboratory & European Spallation Source

This file is part of Mantid.

Mantid is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Mantid is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

File change history is stored at: <https://github.com/mantidproject/mantid>
Code Documentation is available at: <http://doxygen.mantidproject.org>
*/

// Forward deceleration of structs defined at end of header

/** Enum for use in EnabledWhenProperty */
enum ePropertyCriterion {
  IS_DEFAULT,
  IS_NOT_DEFAULT,
  IS_EQUAL_TO,
  IS_NOT_EQUAL_TO,
  IS_MORE_OR_EQ
};

/** Enum for use when combining two EnabledWhenPropertyItems */
enum eLogicOperator { AND, OR, XOR };

class DLLExport EnabledWhenProperty : public IPropertySettings {
public:
  /// Constructs a EnabledWhenProperty object which checks the property
  /// with name given and if it matches the criteria enables it
  EnabledWhenProperty(const std::string &otherPropName,
                      const ePropertyCriterion when,
                      const std::string &value = "");

  /// Constructs a EnabledWhenProperty object which takes ownership of two
  /// already constructed EnabledWhenProperty objects and returns the result
  /// of both of them with the specified logic operator
  EnabledWhenProperty(std::unique_ptr<EnabledWhenProperty> &&conditionOne,
                      std::unique_ptr<EnabledWhenProperty> &&conditionTwo,
                      eLogicOperator logicalOperator);

  /// Copy constructor for EnabledWhenProperty
  EnabledWhenProperty(const EnabledWhenProperty &original);

  /// Checks two EnabledWhenProperty objects match the logic operator
  /// specified and returns the result of both of them
  virtual bool checkComparison(const IPropertyManager *algo) const;

  /// Checks that the specified property matches the criteria given
  virtual bool checkCriterion(const IPropertyManager *algo) const;

  /// Return true/false based on whether the other property satisfies the
  /// criterion
  bool isEnabled(const IPropertyManager *algo) const override;

  /// Return true always
  bool isVisible(const IPropertyManager *algo) const override;

  /// Stub function to satisfy the interface.
  void modify_allowed_values(Property *const);

  /// Make a copy of the present type of validator
  IPropertySettings *clone() override;

protected:
  /// Struct which holds associated property details for comparison
  struct PropertyDetails {
    /// Constructor
    // Have to include this for make_unique to be able to forward arguments
    PropertyDetails(const std::string &otherPropName,
                    const ePropertyCriterion criterion,
                    const std::string &value)
        : otherPropName(otherPropName), criterion(criterion), value(value) {}
    /// Name of the OTHER property that we will check.
    const std::string otherPropName;
    /// Criterion to evaluate
    const ePropertyCriterion criterion;
    /// For the IS_EQUAL_TO or IS_NOT_EQUAL_TO condition,
    /// the value (as string) to check for
    const std::string value;
  };

  /// Struct which holds details for comparison between
  /// two EnabledWhenPropertyObjects
  struct ComparisonDetails {
    /// Constructor
    ComparisonDetails(std::unique_ptr<EnabledWhenProperty> conditionOne,
                      std::unique_ptr<EnabledWhenProperty> conditionTwo,
                      eLogicOperator logicOperator)
        : conditionOne(std::move(conditionOne)),
          conditionTwo(std::move(conditionTwo)), logicOperator(logicOperator) {}

    /// Copy constructor
    ComparisonDetails(const ComparisonDetails &original)
        : conditionOne{new EnabledWhenProperty{*original.conditionOne}},
          conditionTwo{new EnabledWhenProperty{*original.conditionTwo}},
          logicOperator{original.logicOperator} {}

    std::unique_ptr<EnabledWhenProperty> conditionOne;
    std::unique_ptr<EnabledWhenProperty> conditionTwo;
    const eLogicOperator logicOperator;
  };

  /// Checks the algorithm and property are both valid and attempts
  /// to get the value associated with the property
  std::string getPropertyValue(const IPropertyManager *algo) const;

  /// Holds the various details used within the comparison
  std::unique_ptr<PropertyDetails> m_propertyDetails = nullptr;
  /// Holds an object containing details of multiple comparisons
  std::unique_ptr<ComparisonDetails> m_comparisonDetails = nullptr;
};

} // namespace Kernel
} // namespace Mantid

#endif /* MANTID_KERNEL_ENABLEDWHENPROPERTY_H_ */
