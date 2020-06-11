// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidCurveFitting/Algorithms/DoublePulseFit.h"
#include "MantidCurveFitting/CostFunctions/CostFuncFitting.h"

#include "MantidAPI/CompositeFunction.h"
#include "MantidAPI/FuncMinimizerFactory.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/IFuncMinimizer.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/TableRow.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/WorkspaceGroup.h"
#include "MantidCurveFitting/Functions/DeltaFunction.h"

#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/Exception.h"
#include "MantidKernel/StartsWithValidator.h"

#include <math.h>
#include <memory>

using namespace Mantid::CurveFitting::Functions;

namespace Mantid {
namespace CurveFitting {
namespace Algorithms {

void setMultiDataProperties(const Mantid::API::IAlgorithm &qensFit,
                            Mantid::API::IAlgorithm &fit,
                            const Mantid::API::MatrixWorkspace_sptr &workspace,
                            const std::string &suffix) {
  fit.setProperty("InputWorkspace" + suffix, workspace);

  int workspaceIndex = qensFit.getProperty("WorkspaceIndex" + suffix);
  fit.setProperty("WorkspaceIndex" + suffix, workspaceIndex);

  double startX = qensFit.getProperty("StartX" + suffix);
  double endX = qensFit.getProperty("EndX" + suffix);
  fit.setProperty("StartX" + suffix, startX);
  fit.setProperty("EndX" + suffix, endX);

  std::vector<double> exclude = qensFit.getProperty("Exclude" + suffix);
  fit.setProperty("Exclude" + suffix, exclude);
}

void setMultiDataProperties(
    const Mantid::API::IAlgorithm &qensFit, Mantid::API::IAlgorithm &fit,
    const std::vector<Mantid::API::MatrixWorkspace_sptr> &workspaces) {
  setMultiDataProperties(qensFit, fit, workspaces[0], "");

  for (auto i = 1u; i < workspaces.size(); ++i)
    setMultiDataProperties(qensFit, fit, workspaces[i],
                           "_" + std::to_string(i));
}

Mantid::API::IFunction_sptr
getDoublePulseFunction(std::shared_ptr<const API::ParamFunction> const function,
                       double offset, double firstPulseWeight,
                       double secondPulseWeight) {
  auto clonedFunction = function->clone();

  // const double muon_halflife = 2.2;
  // double decay = exp(-offset / muon_halflife);
  // double first_pulse_weighting = decay / (1 + decay);
  // double second_pulse_weighting = 1 / (1 + decay);

  auto convolution =
      std::make_shared<Mantid::CurveFitting::Functions::Convolution>();
  auto delta1 =
      std::make_shared<Mantid::CurveFitting::Functions::DeltaFunction>();
  auto delta2 =
      std::make_shared<Mantid::CurveFitting::Functions::DeltaFunction>();
  auto deltaComposite = std::make_shared<Mantid::API::CompositeFunction>();

  convolution->setAttributeValue("FixResolution", false);
  delta1->setParameter("Centre", -offset);
  delta1->setParameter("Height", firstPulseWeight);
  delta1->fixAll();
  delta2->setParameter("Centre", 0.0);
  delta2->setParameter("Height", secondPulseWeight);
  delta2->fixAll();
  deltaComposite->addFunction(delta1);
  deltaComposite->addFunction(delta2);

  convolution->addFunction(clonedFunction);
  convolution->addFunction(deltaComposite);

  return convolution;
}

Mantid::API::IFunction_sptr getDoublePulseFunction(
    std::shared_ptr<const API::MultiDomainFunction> const initialFunction,
    double offset, double firstPulseWeight, double secondPulseWeight) {
  auto doublePulseFunc = std::make_shared<API::MultiDomainFunction>();
  for (size_t domain = 0; domain < initialFunction->getNumberDomains();
       domain++) {
    auto innerFunction = std::dynamic_pointer_cast<API::ParamFunction>(
        initialFunction->getFunction(domain));
    if (innerFunction) {
      auto twoPulseInnerFunction = getDoublePulseFunction(
          innerFunction, offset, firstPulseWeight, secondPulseWeight);
      doublePulseFunc->addFunction(twoPulseInnerFunction);
      doublePulseFunc->setDomainIndex(domain, domain);
    } else {
      throw std::runtime_error("Convolution does not support composite "
                               "functions. DoublePulseFit.cpp");
    }
  }
  return doublePulseFunc;
}

Mantid::API::IFunction_sptr extractInnerFunction(
    std::shared_ptr<const Mantid::CurveFitting::Functions::Convolution> const
        convolutionFunction) {
  return convolutionFunction->getFunction(0);
}

Mantid::API::IFunction_sptr extractInnerFunction(
    std::shared_ptr<const API::MultiDomainFunction> const doublePulseFunc) {
  auto extractedFunction = std::make_shared<API::MultiDomainFunction>();
  for (size_t domain = 0; domain < doublePulseFunc->getNumberDomains();
       domain++) {
    auto convFunction = std::dynamic_pointer_cast<Convolution>(
        doublePulseFunc->getFunction(domain));
    if (convFunction) {
      extractedFunction->addFunction(extractInnerFunction(convFunction));
      extractedFunction->setDomainIndex(domain, domain);
    } else {
      throw std::runtime_error(
          "Cannot extract from non convolution function. DoublePulseFit.cpp");
    }
  }
  return extractedFunction;
}

// Register the class into the algorithm factory
DECLARE_ALGORITHM(DoublePulseFit)

/// Default constructor
DoublePulseFit::DoublePulseFit() : IFittingAlgorithm() {}

/** Initialisation method
 */
void DoublePulseFit::initConcrete() {
  declareProperty("Ties", "", Kernel::Direction::Input);
  getPointerToProperty("Ties")->setDocumentation(
      "Math expressions defining ties between parameters of "
      "the fitting function.");
  declareProperty("Constraints", "", Kernel::Direction::Input);
  getPointerToProperty("Constraints")->setDocumentation("List of constraints");
  auto mustBePositive = std::make_shared<Kernel::BoundedValidator<int>>();
  mustBePositive->setLower(0);
  declareProperty(
      "MaxIterations", 500, mustBePositive->clone(),
      "Stop after this number of iterations if a good fit is not found");
  declareProperty("OutputStatus", "", Kernel::Direction::Output);
  getPointerToProperty("OutputStatus")
      ->setDocumentation("Whether the fit was successful");
  declareProperty("OutputChi2overDoF", 0.0, "Returns the goodness of the fit",
                  Kernel::Direction::Output);

  std::vector<std::string> minimizerOptions =
      API::FuncMinimizerFactory::Instance().getKeys();
  Kernel::IValidator_sptr minimizerValidator =
      std::make_shared<Kernel::StartsWithValidator>(minimizerOptions);

  declareProperty("Minimizer", "Levenberg-Marquardt", minimizerValidator,
                  "Minimizer to use for fitting.");

  std::vector<std::string> costFuncOptions =
      API::CostFunctionFactory::Instance().getKeys();
  // select only CostFuncFitting variety
  for (auto &costFuncOption : costFuncOptions) {
    auto costFunc = std::dynamic_pointer_cast<CostFunctions::CostFuncFitting>(
        API::CostFunctionFactory::Instance().create(costFuncOption));
    if (!costFunc) {
      costFuncOption = "";
    }
  }
  Kernel::IValidator_sptr costFuncValidator =
      std::make_shared<Kernel::ListValidator<std::string>>(costFuncOptions);
  declareProperty(
      "CostFunction", "Least squares", costFuncValidator,
      "The cost function to be used for the fit, default is Least squares",
      Kernel::Direction::InOut);
  declareProperty(
      "CreateOutput", false,
      "Set to true to create output workspaces with the results of the fit"
      "(default is false).");
  declareProperty(
      "Output", "",
      "A base name for the output workspaces (if not "
      "given default names will be created). The "
      "default is to use the name of the original data workspace as prefix "
      "followed by suffixes _Workspace, _Parameters, etc.");
  declareProperty("CalcErrors", false,
                  "Set to true to calcuate errors when output isn't created "
                  "(default is false).");
  declareProperty("OutputCompositeMembers", false,
                  "If true and CreateOutput is true then the value of each "
                  "member of a Composite Function is also output.");
  declareProperty(std::make_unique<Kernel::PropertyWithValue<bool>>(
                      "ConvolveMembers", false),
                  "If true and OutputCompositeMembers is true members of any "
                  "Convolution are output convolved\n"
                  "with corresponding resolution");
  declareProperty("OutputParametersOnly", false,
                  "Set to true to output only the parameters and not "
                  "workspace(s) with the calculated values\n"
                  "(default is false, ignored if CreateOutput is false and "
                  "Output is an empty string).");
  declareProperty("PulseOffset", 0.0, "The time offset between the two pulses");
  declareProperty("FirstPulseWeight", 1.0, "Weighting of first pulse.");
  declareProperty("SecondPulseWeight", 1.0, "Weighting of first pulse.");
}

void DoublePulseFit::execConcrete() {
  setOutputProperties();
  declareAdditionalProperties();

  Mantid::API::Algorithm_sptr fitAlg = createChildAlgorithm("Fit");
  Mantid::API::IFunction_sptr function = getProperty("Function");
  Mantid::API::IFunction_sptr doublePulseFunction;
  auto pulseOffset = getProperty("PulseOffset");
  auto firstPulseWeight = getProperty("FirstPulseWeight");
  auto secondPulseWeight = getProperty("SecondPulseWeight");
  if (auto paramFunction =
          std::dynamic_pointer_cast<Mantid::API::ParamFunction>(function)) {
    doublePulseFunction = getDoublePulseFunction(
        paramFunction, pulseOffset, firstPulseWeight, secondPulseWeight);
  } else if (auto multiDomainFunction =
                 std::dynamic_pointer_cast<Mantid::API::MultiDomainFunction>(
                     function)) {
    doublePulseFunction = getDoublePulseFunction(
        multiDomainFunction, pulseOffset, firstPulseWeight, secondPulseWeight);
  } else {
    throw std::runtime_error("Incompatible function form. DoublePulseFit.cpp");
  }

  runFitAlgorith(fitAlg, doublePulseFunction, getProperty("MaxIterations"));

  createOutput(fitAlg, doublePulseFunction);
}

/**
 * This function declares the additional output properties which are not
 * present in init. It is done this way to match the behaviour of Fit which
 * declares many of its output properties internally setting the names of
 * output workspaces based upon the input property Output.
 */
void DoublePulseFit::declareAdditionalProperties() {
  std::string baseName = getProperty("Output");
  if (baseName.empty()) {
    API::Workspace_const_sptr ws = getProperty("InputWorkspace");
    baseName = ws->getName();
  }

  Mantid::API::IFunction_sptr function = getProperty("Function");
  if (m_makeOutput) {
    declareProperty(
        std::make_unique<API::WorkspaceProperty<API::ITableWorkspace>>(
            "OutputNormalisedCovarianceMatrix", "", Kernel::Direction::Output),
        "The name of the TableWorkspace in which to store the final "
        "covariance "
        "matrix");
    setPropertyValue("OutputNormalisedCovarianceMatrix",
                     baseName + "NormalisedCovarianceMatrix");

    declareProperty(
        std::make_unique<API::WorkspaceProperty<API::ITableWorkspace>>(
            "OutputParameters", "", Kernel::Direction::Output),
        "The name of the TableWorkspace in which to store the "
        "final fit parameters");

    setPropertyValue("OutputParameters", baseName + "Parameters");

    if (m_outputFitData) {
      if (m_multiDomain) {
        declareProperty(
            std::make_unique<
                API::WorkspaceProperty<Mantid::API::WorkspaceGroup>>(
                "OutputWorkspace", "", Kernel::Direction::Output),
            "Name of the output Workspace holding resulting simulated "
            "spectrum");
        setPropertyValue("OutputWorkspace", baseName + "Workspace");
      } else {
        declareProperty(
            std::make_unique<
                API::WorkspaceProperty<Mantid::API::MatrixWorkspace>>(
                "OutputWorkspace", "", Kernel::Direction::Output),
            "Name of the output Workspace holding resulting simulated "
            "spectrum");
        setPropertyValue("OutputWorkspace", baseName + "Workspace");
      }
    }
  }
}

std::vector<Mantid::API::MatrixWorkspace_sptr>
DoublePulseFit::getWorkspaces() const {
  std::vector<Mantid::API::MatrixWorkspace_sptr> workspaces;
  workspaces.reserve(m_workspacePropertyNames.size());
  for (const auto &propertyName : m_workspacePropertyNames) {
    Mantid::API::Workspace_sptr workspace = getProperty(propertyName);
    workspaces.emplace_back(
        std::dynamic_pointer_cast<Mantid::API::MatrixWorkspace>(workspace));
  }
  return workspaces;
}

void DoublePulseFit::setOutputProperties() {
  std::string baseName = getProperty("Output");
  bool makeOutput = getProperty("CreateOutput");
  m_makeOutput = makeOutput || !baseName.empty();
  m_outputFitData = !getProperty("OutputParametersOnly");
  Mantid::API::IFunction_sptr function = getProperty("Function");
  m_multiDomain = function->getNumberDomains() > 1;
}

void DoublePulseFit::runFitAlgorith(Mantid::API::IAlgorithm_sptr fitAlg,
                                    Mantid::API::IFunction_sptr function,
                                    int maxIterations) {
  const bool convolveMembers = getProperty("ConvolveMembers");
  const bool ignoreInvalidData = getProperty("IgnoreInvalidData");
  const bool calcErrors = getProperty("CalcErrors");
  const bool createOutput = getProperty("CreateOutput");
  const bool outputCompositeMembers = getProperty("OutputCompositeMembers");
  const bool outputParametersOnly = getProperty("OutputParametersOnly");

  fitAlg->setProperty("Function", function);
  setMultiDataProperties(*this, *fitAlg, getWorkspaces());
  fitAlg->setProperty("IgnoreInvalidData", ignoreInvalidData);
  fitAlg->setProperty("DomainType", getPropertyValue("DomainType"));
  fitAlg->setProperty("EvaluationType", getPropertyValue("EvaluationType"));
  fitAlg->setPropertyValue("PeakRadius", getPropertyValue("PeakRadius"));
  fitAlg->setProperty("Ties", getPropertyValue("Ties"));
  fitAlg->setProperty("Constraints", getPropertyValue("Constraints"));
  fitAlg->setProperty("MaxIterations", maxIterations);
  fitAlg->setProperty("Minimizer", getPropertyValue("Minimizer"));
  fitAlg->setProperty("CostFunction", getPropertyValue("CostFunction"));
  fitAlg->setProperty("CalcErrors", calcErrors);
  fitAlg->setProperty("OutputCompositeMembers", outputCompositeMembers);
  fitAlg->setProperty("ConvolveMembers", convolveMembers);
  fitAlg->setProperty("CreateOutput", createOutput);
  fitAlg->setProperty("OutputParametersOnly", outputParametersOnly);
  fitAlg->setProperty("Output", getPropertyValue("Output"));
  fitAlg->executeAsChildAlg();
}

void DoublePulseFit::createOutput(Mantid::API::IAlgorithm_sptr fitAlg,
                                  Mantid::API::IFunction_sptr function) {
  Mantid::API::IFunction_sptr extractedFunction;
  if (auto convFunction = std::dynamic_pointer_cast<
          Mantid::CurveFitting::Functions::Convolution>(function)) {
    extractedFunction = extractInnerFunction(convFunction);
  } else if (auto multiDomainFunction =
                 std::dynamic_pointer_cast<Mantid::API::MultiDomainFunction>(
                     function)) {
    extractedFunction = extractInnerFunction(multiDomainFunction);
  } else {
    throw std::runtime_error("Incompatible function form. DoublePulseFit.cpp");
  }
  double outputChi2overDoF = fitAlg->getProperty("OutputChi2overDoF");
  std::string outputStatus = fitAlg->getProperty("OutputStatus");
  setProperty("OutputStatus", outputStatus);
  setProperty("OutputChi2overDoF", outputChi2overDoF);
  if (m_makeOutput) {
    // The scientists want the details of the double pulse analysis to be
    // hidden.
    // This alg is used to produce the output covariance matrix and property
    // table as though the double pulse functions were not there.
    Mantid::API::Algorithm_sptr figAlgForOutputs = createChildAlgorithm("Fit");
    runFitAlgorith(figAlgForOutputs, extractedFunction, 0);

    Mantid::API::ITableWorkspace_sptr covarianceMatrix =
        figAlgForOutputs->getProperty("OutputNormalisedCovarianceMatrix");
    setProperty("OutputNormalisedCovarianceMatrix", covarianceMatrix);
    Mantid::API::ITableWorkspace_sptr parametersTable =
        figAlgForOutputs->getProperty("OutputParameters");
    setProperty("OutputParameters", parametersTable);
    if (m_outputFitData) {
      if (m_multiDomain) {
        Mantid::API::WorkspaceGroup_sptr outputWorkspace =
            fitAlg->getProperty("OutputWorkspace");
        setProperty("OutputWorkspace", outputWorkspace);
      } else {
        Mantid::API::MatrixWorkspace_sptr outputWorkspace =
            fitAlg->getProperty("OutputWorkspace");
        setProperty("OutputWorkspace", outputWorkspace);
      }
    }
  }
  //   // Fit has the function property as an InOut object and modifies it in
  //   place.
  //   // To replicate this behaviour we reset the relevant pointer to the new
  //   // function.
  //   Mantid::API::IFunction_sptr inputFunction = getProperty("Function");
  //   inputFunction.reset(extractedFunction.get());
  //   extractedFunction.reset();
  setProperty("Function", extractedFunction);
}

} // namespace Algorithms
} // namespace CurveFitting
} // namespace Mantid