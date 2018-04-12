#ifndef MANTID_ALGORITHMS_QENSFITSEQUENTIAL_H_
#define MANTID_ALGORITHMS_QENSFITSEQUENTIAL_H_

#include "MantidAPI/Algorithm.h"
#include "MantidAPI/ITableWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceGroup.h"
#include "MantidKernel/IValidator.h"

namespace Mantid {
namespace Algorithms {

/**
  Copyright &copy; 2015 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
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
class DLLExport QENSFitSequential : public API::Algorithm {
public:
  const std::string name() const override;
  int version() const override;
  const std::string category() const override;
  const std::string summary() const override;

protected:
  virtual std::vector<API::MatrixWorkspace_sptr> getWorkspaces() const;
  virtual void deleteTemporaryWorkspaces(const std::string &outputBaseName);
  virtual API::ITableWorkspace_sptr performFit(const std::string &input);

private:
  void init() override;
  virtual void initConcrete();
  virtual void setup();
  void exec() override;
  virtual void postExec(MatrixWorkspace_sptr result);
  virtual Kernel::IValidator_sptr functionValidator() const;

  std::string getInputString(
      const std::vector<API::MatrixWorkspace_sptr> &workspaces) const;
  virtual std::vector<std::string> getFitParameterNames() const;
  API::MatrixWorkspace_sptr
  processIndirectFitParameters(API::ITableWorkspace_sptr parameterWorkspace);

  void renameWorkspaces(API::WorkspaceGroup_sptr outputGroup,
                        const std::vector<std::string> &spectra);
  API::MatrixWorkspace_sptr
  copyLogs(API::MatrixWorkspace_sptr resultWorkspace,
           const std::vector<API::MatrixWorkspace_sptr> &workspaces);
  void extractMembers(API::WorkspaceGroup_sptr resultGroupWs,
                      const std::vector<API::MatrixWorkspace_sptr> &workspaces,
                      const std::string &outputWsName);

  API::IAlgorithm_sptr
  extractMembersAlgorithm(API::WorkspaceGroup_sptr resultGroupWs,
                          const std::string &outputWsName) const;
};

} // namespace Algorithms
} // namespace Mantid

#endif
