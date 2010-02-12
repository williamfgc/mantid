#ifndef MANTID_ALGORITHMS_FLATBACKGROUND_H_
#define MANTID_ALGORITHMS_FLATBACKGROUND_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/Algorithm.h"

namespace Mantid
{
namespace Algorithms
{
/** Finds a constant value fit to an appropriate range of each desired spectrum
    and subtracts that value from the entire spectrum.

    Required Properties:
    <UL>
    <LI> InputWorkspace     - The name of the input workspace. </LI>
    <LI> OutputWorkspace    - The name to give the output workspace. </LI>
    <LI> SpectrumIndexList  - The workspace indices of the spectra to fit background to. </LI>
    <LI> StartX             - The start of the flat region to fit to. </LI>
    <LI> EndX               - The end of the flat region to fit to. </LI>
    <LI> Mode               - How to estimate the background number of counts: a linear fit or the mean. </LI>
    </UL>

    @author Russell Taylor, Tessella plc
    @date 5/02/2009

    Copyright &copy; 2009 STFC Rutherford Appleton Laboratory

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

    File change history is stored at: <https://svn.mantidproject.org/mantid/trunk/Code/Mantid>
    Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class DLLExport FlatBackground : public API::Algorithm
{
public:
  /// (Empty) Constructor
  FlatBackground() : API::Algorithm(), m_progress(NULL) {}
  /// Virtual destructor
  virtual ~FlatBackground() {if(m_progress) delete m_progress;m_progress=NULL;}
  /// Algorithm's name
  virtual const std::string name() const { return "FlatBackground"; }
  /// Algorithm's version
  virtual const int version() const { return (1); }
  /// Algorithm's category for identification
  virtual const std::string category() const { return "SANS"; }

private:
  /// Initialisation code
  void init();
  ///Execution code
  void exec();

  void checkRange(double& startX, double& endX);
  bool isModeMean(const std::string &mode);
  void getSpecInds(std::vector<int> &output, const int workspaceTotal);
  double Mean(const API::MatrixWorkspace_const_sptr WS, const int specInd, const double startX, const double endX, double &variance) const;
  double LinearFit(API::MatrixWorkspace_sptr WS, int spectrum, double startX, double endX);

  /// Progress reporting
  API::Progress* m_progress;

};

} // namespace Algorithms
} // namespace Mantid

#endif /*MANTID_ALGORITHMS_FLATBACKGROUND_H_*/
