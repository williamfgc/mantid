#ifndef MANTIDQTCUSTOMINTERFACES_MUONANALYSIS_H_
#define MANTIDQTCUSTOMINTERFACES_MUONANALYSIS_H_

//----------------------
// Includes
//----------------------
#include "ui_MuonAnalysis.h"
#include "MantidQtAPI/UserSubWindow.h"

#include "MantidQtMantidWidgets/pythonCalc.h"
#include "MantidQtMantidWidgets/MWDiag.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/MatrixWorkspace.h"
#include <map>

namespace MantidQt
{
namespace CustomInterfaces
{

namespace Muon
{
  class MuonAnalysisOptionTab;
  class MuonAnalysisFitDataTab;
}


/** 
This is the main class for the MuonAnalysis interface
see <http://www.mantidproject.org/MuonAnalysis>.    

@author Anders Markvardsen, ISIS, RAL

Copyright &copy; 2010 ISIS Rutherford Appleton Laboratory & NScD Oak Ridge National Laboratory

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


class MuonAnalysis : public MantidQt::API::UserSubWindow
{
  Q_OBJECT

public:
  /// Name of the interface
  static std::string name() { return "Muon Analysis"; }

public:
  /// Default Constructor
  MuonAnalysis(QWidget *parent = 0);

private slots:
  /// Exit the interface
  //void exitClicked();

  /// Guess Alpha clicked
  void guessAlphaClicked();

  /// Input file changed in MWRunFiles widget
  void inputFileChanged_MWRunFiles();

  // Load current button
  void runLoadCurrent();

  /// group table changed
  void groupTableChanged(int row, int column);

  // group table clicked
  void groupTableClicked(int row, int column);

  // group table vertical label clicked
  void groupTableClicked(int row);

  /// group table changed
  void pairTableChanged(int row, int column);

  // pair table clicked
  void pairTableClicked(int row, int column);

  // pair table vertical lable clicked
  void pairTableClicked(int row);

  /// group table plot button
  void runGroupTablePlotButton();

  /// group table plot button
  void runPairTablePlotButton();

  /// Save grouping button
  void runSaveGroupButton();

  /// Load grouping button
  void runLoadGroupButton();

  /// Clear grouping button
  void runClearGroupingButton(); 

  /// User select instrument
  void userSelectInstrument(const QString& prefix);

  ///
  void runFrontPlotButton();

  ///
  void runFrontGroupGroupPairComboBox(int index);

  ///
  void muonAnalysisHelpClicked();

  ///
  void muonAnalysisHelpGroupingClicked();

  ///
  void muonAnalysisHelpPlottingClicked();

  ///
  void muonAnalysisHelpDataAnalysisClicked();
  
  ///
  void runFirstGoodBinFront();

  /// Check to see if the user want to append the previous run and set accordingly
  void checkAppendingPreviousRun();

  /// Check to see if the user want to append the next run and set accordingly
  void checkAppendingNextRun();

  /// Changes the format of the file string. i.e file MUSR001-06 are files 01-06
  ///void appendSelected();

  /// Bunch the data (rebin)
  void reBunch(const std::string&);

  /// Remove any rebinned data so that the plot shows the raw data again
  void makeRaw(const std::string&);

  ///
  void changeTab(int);

  /// Assigns a peak picker tool to the workspace (@param::workspace name)
  void assignPeakPickerTool(const QString &);

  /// Change the fit style and color
  void changeFitPlotType(const QString &);

  /// Change the data style and color
  void changeDataPlotType(const QString &);

  void groupFittedWorkspaces(QString);

  /// This slot is a placeholder for executing MuonAnalysisFitDataTab::beforeDoFit 
  void beforeDoFit(const QtBoolPropertyManager*);

private:
  /// Initialize the layout
  virtual void initLayout();

  /// Set start up interface look
  void startUpLook();

  /// Input file changed - update GUI accordingly
  void inputFileChanged(const QString& filename);

  /// Return the pair which is in focus and -1 if none
  int pairInFocus();

  /// is grouping set
  bool isGroupingSet();

  /// Apply grouping specified in xml file to workspace
  bool applyGroupingToWS( const std::string& inputWS,  const std::string& outputWS, 
    const std::string& filename);

  /// create WS contained the data for a plot
  void createPlotWS(const std::string& groupName, const std::string& wsname);

  /// Apply whatever grouping is specified in GUI tables to workspace
  bool applyGroupingToWS( const std::string& inputWS,  const std::string& outputWS);

  /// Update front 
  void updateFront();

  /// Update front anc pair combo box
  void updateFrontAndCombo();

  /// Calculate number of detectors from string of type 1-3, 5, 10-15
  int numOfDetectors(const std::string& str) const;

  /// Return a vector of IDs for row number from string of type 1-3, 5, 10-15
  std::vector<int> spectrumIDs(const std::string& str) const;

  /// is string a number?
  bool isNumber(const std::string& s) const;

  /// Clear tables and front combo box
  void clearTablesAndCombo();

  /// When no data loaded set various buttons etc to inactive
  void noDataAvailable();

  /// When data loaded set various buttons etc to active
  void nowDataAvailable();

  /// Check if grouping in table is consistent with data file
  std::string isGroupingAndDataConsistent();

  ///Return true if data are loaded
  bool areDataLoaded();

  /// Return number of pairs
  int numPairs();

  /// Return number of groups defined (not including pairs)
  int numGroups();

  /// Plot group
  void plotGroup(const std::string& plotType);

  /// Plot pair
  void plotPair(const std::string& plotType);

  //The form generated by Qt Designer
  Ui::MuonAnalysis m_uiForm;

  /// group plot functions
  QStringList m_groupPlotFunc;

  /// pair plot functions
  QStringList m_pairPlotFunc;

  /// The last directory that was viewed
  QString m_last_dir;

  /// name of workspace
  std::string m_workspace_name;

  /// name of the loaded data
  QString m_currentDataName;

  /// boolean to tell whether the fit property browser has been assigned
  bool m_assigned;

  /// which group table row has the user last clicked on
  int m_groupTableRowInFocus;

  /// which pair table row has the user last clicked on
  int m_pairTableRowInFocus;

  /// Index of the current tab.
  int m_tabNumber;

  /// used to test that a new filename has been entered 
  QString m_previousFilename;

  /// used to test that a new number of steps to rebin by has been entered
  QString m_previousRebinSteps;

  /// used to test whether the data has already been bunched before
  std::string m_previousBunchWsName;

  /// List of current group names 
  std::vector<std::string> m_groupNames;

  /// name for file to temperary store grouping
  std::string m_groupingTempFilename;

  ///
  void updatePairTable();

  /// Currently selected instrument
  QString m_curInterfaceSetup;

  /// tell which group is in which row
  std::vector<int> m_pairToRow;

  /// tell which group is in which row
  std::vector<int> m_groupToRow;

  ///
  void checkIf_ID_dublicatesInTable(const int row);

  /// Return the group-number for the group in a row. 
  /// Return -1 if invalid group in row
  int getGroupNumberFromRow(int row);

  /// Return the pair-number for the pair in a row. 
  /// Return -1 if invalid pair in row
  int getPairNumberFromRow(int row);

  /// first good bin returned in ms
  QString firstGoodBin();

  /// According to Plot Options what is the time to plot from in ms
  double plotFromTime();

  /// According to Plot Options what is the time to plot to in ms
  double plotToTime();

  /// time zero returned in ms
  QString timeZero();

  /// set grouping in table from information from nexus raw file
  void setGroupingFromNexus(const QString& nexusFile); 

  ///
  void setDummyGrouping(const int numDetectors);

  ///
  void setGroupingFromIDF(const std::string& mainFieldDirection, Mantid::API::MatrixWorkspace_sptr matrix_workspace);

  /// title of run
  std::string m_title;

  /// group defaults are saved to
  QString m_settingsGroup;

  /// Load auto saved values
  void loadAutoSavedValues(const QString& group);

  /// connect the settings for the fit function to their respective slots
  void loadFittings();

  /// Add or take one away from the range of files
  void setAppendingRun(int inc);

  /// change and load the run depending on the value passed as a parameter
  void changeRun(int amountToChange);

  /// Separate the muon file, current File becomes the beginning part of the file i.e the instrument in the file MUSR002413
  void separateMuonFile(QString & currentFile, int & runNumberSize, int & runNumber);

  /// Include the 0's fromt eh beginning of the file that were lost in conversion from QString to int
  void getFullCode(int size, QString & limitedCode);

  /// handles option tab work
  MantidQt::CustomInterfaces::Muon::MuonAnalysisOptionTab* m_optionTab;
  /// handles fit data work
  MantidQt::CustomInterfaces::Muon::MuonAnalysisFitDataTab* m_fitDataTab;

  //A reference to a logger
  static Mantid::Kernel::Logger & g_log;
};

}
}

#endif //MANTIDQTCUSTOMINTERFACES_MUONANALYSIS_H_
