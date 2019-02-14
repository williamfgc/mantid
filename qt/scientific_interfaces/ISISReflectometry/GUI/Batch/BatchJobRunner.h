// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_CUSTOMINTERFACES_BATCHJOBRUNNER_H_
#define MANTID_CUSTOMINTERFACES_BATCHJOBRUNNER_H_

#include "Common/DllConfig.h"
#include "IBatchJobRunner.h"
#include "MantidAPI/IAlgorithm_fwd.h"
#include "MantidQtWidgets/Common/BatchAlgorithmRunner.h"
#include "Reduction/Batch.h"

namespace MantidQt {
namespace CustomInterfaces {

// Override the configured algorithm to pass to BatchAlgorithmRunner
// so that we can associate our own data with it
class BatchJobAlgorithm : public MantidQt::API::ConfiguredAlgorithm {
public:
  BatchJobAlgorithm(
      Mantid::API::IAlgorithm_sptr algorithm,
      MantidQt::API::ConfiguredAlgorithm::AlgorithmRuntimeProps properties,
      Item *item);

  Item *item();

private:
  // The data is an item in the table (i.e. a row or group)
  Item *m_item;
};

using BatchJobAlgorithm_sptr = boost::shared_ptr<BatchJobAlgorithm>;

/**
 * The BatchJobRunner class manages the running of algorithms based
 * on the reduction configuration and handling of state when algorithms
 * complete
 */
class MANTIDQT_ISISREFLECTOMETRY_DLL BatchJobRunner : public IBatchJobRunner {
public:
  explicit BatchJobRunner(Batch batch);

  bool isProcessing() const override;
  bool isAutoreducing() const override;

  void resumeReduction() override;
  void reductionPaused() override;
  void resumeAutoreduction() override;
  void autoreductionPaused() override;

  void setReprocessFailedItems(bool reprocessFailed) override;

  void
  algorithmStarted(MantidQt::API::ConfiguredAlgorithm_sptr algorithm) override;
  void
  algorithmComplete(MantidQt::API::ConfiguredAlgorithm_sptr algorithm) override;
  void algorithmError(MantidQt::API::ConfiguredAlgorithm_sptr algorithm,
                      std::string const &message) override;

  std::vector<std::string> algorithmOutputWorkspacesToSave(
      MantidQt::API::ConfiguredAlgorithm_sptr algorithm) const override;

  void notifyWorkspaceDeleted(std::string const &wsName) override;
  void notifyWorkspaceRenamed(std::string const &oldName,
                              std::string const &newName) override;
  void notifyAllWorkspacesDeleted() override;

  std::deque<MantidQt::API::ConfiguredAlgorithm_sptr> getAlgorithms() override;

private:
  Batch m_batch;
  bool m_isProcessing;
  bool m_isAutoreducing;
  bool m_reprocessFailed;
  bool m_processAll;

  std::vector<std::string> getWorkspacesToSave(Group const &group) const;
  std::vector<std::string> getWorkspacesToSave(Row const &row) const;
};
} // namespace CustomInterfaces
} // namespace MantidQt

#endif // MANTID_CUSTOMINTERFACES_BATCHJOBRUNNER_H_
