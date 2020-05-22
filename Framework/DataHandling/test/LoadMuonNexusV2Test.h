// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2020 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/Run.h"
#include "MantidAPI/WorkspaceGroup.h"
#include "MantidDataHandling/LoadMuonNexusV2.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidKernel/UnitLabelTypes.h"

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <vector>

using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::DataHandling;
using namespace Mantid::DataObjects;

class LoadMuonNexusV2Test : public CxxTest::TestSuite {
private:
  // helper methods
  std::string createSpectraList(const std::vector<int> &spectraList) {
    std::ostringstream oss;
    std::copy(std::begin(spectraList), std::end(spectraList),
              std::ostream_iterator<int>(oss, ","));
    std::string spectraStringList(oss.str());
    return spectraStringList;
  }

public:
  void tearDown() override { AnalysisDataService::Instance().clear(); }

  void testExec() {
    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");

    TS_ASSERT_THROWS_NOTHING(ld.execute());
    TS_ASSERT(ld.isExecuted());
    // Verify that the output workspace exists
    MatrixWorkspace_sptr output_ws;
    TS_ASSERT_THROWS_NOTHING(
        output_ws = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(
            "outWS"));

    Workspace2D_sptr output2D =
        std::dynamic_pointer_cast<Workspace2D>(output_ws);

    const Mantid::API::Run &run = output2D->run();
    int goodfrm = run.getPropertyAsIntegerValue("goodfrm");
    TS_ASSERT_EQUALS(goodfrm, 14320);
    double firstGoodData = ld.getProperty("FirstGoodData");
    TS_ASSERT_EQUALS(firstGoodData, 0.384);
    double timeZero = ld.getProperty("TimeZero");
    TS_ASSERT_DELTA(timeZero, 0.1599999, 1e-5);

    // Check that timeZero has been applied to the output spectra
    // as LoadISISNexus does not do this.
    // First time reading should be shifted by time zero.
    TS_ASSERT_DELTA(output2D->x(3)[0], -0.1599999, 1e-5);
    TS_ASSERT_DELTA(output2D->x(67)[0], -0.1599999, 1e-5);
    TS_ASSERT_DELTA(output2D->x(81)[0], -0.1599999, 1e-5);

    // Check the unit has been set correctly
    TS_ASSERT_EQUALS(output2D->getAxis(0)->unit()->unitID(), "Label");
    TS_ASSERT(!output2D->isDistribution());

    // Check that sample temp and field set
    double temperature = run.getPropertyAsSingleValue("sample_temp");
    TS_ASSERT_EQUALS(10.0, temperature);
    double field = run.getPropertyAsSingleValue("sample_magn_field");
    TS_ASSERT_EQUALS(20.0, field);
  }

  void testExecWithDeadtimeTable() {

    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");
    const std::string deadTimeWSName = "LoadMuonNexusV2Test_DeadTimes";
    ld.setPropertyValue("DeadTimeTable", deadTimeWSName);

    TS_ASSERT_THROWS_NOTHING(ld.execute());
    TS_ASSERT(ld.isExecuted());
    // Verify that the output workspace exists
    MatrixWorkspace_sptr output_ws;
    TS_ASSERT_THROWS_NOTHING(
        output_ws = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(
            "outWS"));

    Workspace2D_sptr output2D =
        std::dynamic_pointer_cast<Workspace2D>(output_ws);

    // Check detector grouping table
    TableWorkspace_sptr deadTimeTable;
    TS_ASSERT_THROWS_NOTHING(
        deadTimeTable =
            AnalysisDataService::Instance().retrieveWS<TableWorkspace>(
                deadTimeWSName));
    TS_ASSERT(deadTimeTable);
    // Check number of rows and columns
    TS_ASSERT_EQUALS(deadTimeTable->columnCount(), 2);
    TS_ASSERT_EQUALS(deadTimeTable->rowCount(),
                     output2D->getNumberHistograms());
    // Check Deadtimes
    TS_ASSERT_DELTA(deadTimeTable->Double(0, 1), -0.0095861498, 1e-6);
    TS_ASSERT_DELTA(deadTimeTable->Double(20, 1), 0.0067306999, 1e-6);
    TS_ASSERT_DELTA(deadTimeTable->Double(62, 1), 0.0073113599, 1e-6);
  }

  void testExecWithGroupingTable() {

    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");
    const std::string groupingWSName = "LoadMuonNexusV2Test_Grouping";
    ld.setPropertyValue("DetectorGroupingTable", groupingWSName);
    TS_ASSERT_THROWS_NOTHING(ld.execute());
    TS_ASSERT(ld.isExecuted());
    // Verify that the output workspace exists
    MatrixWorkspace_sptr output_ws;
    TS_ASSERT_THROWS_NOTHING(
        output_ws = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(
            "outWS"));
    Workspace2D_sptr output2D =
        std::dynamic_pointer_cast<Workspace2D>(output_ws);

    // Check detector grouping table
    TableWorkspace_sptr groupingTable;
    TS_ASSERT_THROWS_NOTHING(
        groupingTable =
            AnalysisDataService::Instance().retrieveWS<TableWorkspace>(
                groupingWSName));
    TS_ASSERT(groupingTable);
    // Check number of rows and columns
    TS_ASSERT_EQUALS(groupingTable->columnCount(), 1);
    TS_ASSERT_EQUALS(groupingTable->rowCount(), 2);
    // Check grouping
    std::vector<int> testGroupingVec;
    // Half the detectors are in the first group
    for (int i = 1; i < 49; ++i)
      testGroupingVec.emplace_back(i);
    TS_ASSERT_EQUALS(groupingTable->cell<std::vector<int>>(0, 0),
                     testGroupingVec);
    testGroupingVec.clear();
    for (int i = 49; i < static_cast<int>(output2D->getNumberHistograms() + 1);
         ++i)
      testGroupingVec.emplace_back(i);
    TS_ASSERT_EQUALS(groupingTable->cell<std::vector<int>>(1, 0),
                     testGroupingVec);
  }

  void testExecWithSpectraList() {
    std::vector<int> spectraIntegerList = {1, 21, 63};
    auto spectraList = createSpectraList(spectraIntegerList);

    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");
    ld.setPropertyValue("SpectrumList", spectraList);
    const std::string deadTimeWSName = "LoadMuonNexusV2Test_DeadTimes";
    ld.setPropertyValue("DeadTimeTable", deadTimeWSName);
    const std::string groupingWSName = "LoadMuonNexusV2Test_Grouping";
    ld.setPropertyValue("DetectorGroupingTable", groupingWSName);

    TS_ASSERT_THROWS_NOTHING(ld.execute());
    TS_ASSERT(ld.isExecuted());
    // Verify that the output workspace exists
    MatrixWorkspace_sptr output_ws;
    TS_ASSERT_THROWS_NOTHING(
        output_ws = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(
            "outWS"));
    Workspace2D_sptr output2D =
        std::dynamic_pointer_cast<Workspace2D>(output_ws);

    // Test correct spectra loaded
    TS_ASSERT_EQUALS(output2D->getNumberHistograms(), 3);
    // Check that spectra maps to correct detector
    for (size_t i = 0; i < 3; ++i) {
      auto detectorgroup = output2D->getSpectrum(i).getDetectorIDs();
      TS_ASSERT_EQUALS(detectorgroup.size(), 1);
      TS_ASSERT_EQUALS(*detectorgroup.begin(), spectraIntegerList[i]);
    }

    // Check detector grouping table
    TableWorkspace_sptr deadTimeTable;
    TS_ASSERT_THROWS_NOTHING(
        deadTimeTable =
            AnalysisDataService::Instance().retrieveWS<TableWorkspace>(
                deadTimeWSName));
    TS_ASSERT(deadTimeTable);
    // Check number of rows and columns
    TS_ASSERT_EQUALS(deadTimeTable->columnCount(), 2);
    TS_ASSERT_EQUALS(deadTimeTable->rowCount(), 3);
    // Check Deadtimes
    TS_ASSERT_DELTA(deadTimeTable->Double(0, 1), -0.0095861498, 1e-6);
    TS_ASSERT_DELTA(deadTimeTable->Double(1, 1), 0.0067306999, 1e-6);
    TS_ASSERT_DELTA(deadTimeTable->Double(2, 1), 0.0073113599, 1e-6);
  }
  void testExecWithSpectraMax() {

    int specMax = 24;
    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");
    ld.setPropertyValue("SpectrumMax", "24");
    const std::string deadTimeWSName = "LoadMuonNexusV2Test_DeadTimes";
    ld.setPropertyValue("DeadTimeTable", deadTimeWSName);
    const std::string groupingWSName = "LoadMuonNexusV2Test_Grouping";
    ld.setPropertyValue("DetectorGroupingTable", groupingWSName);

    TS_ASSERT_THROWS_NOTHING(ld.execute());
    TS_ASSERT(ld.isExecuted());
    // Verify that the output workspace exists
    MatrixWorkspace_sptr output_ws;
    TS_ASSERT_THROWS_NOTHING(
        output_ws = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>(
            "outWS"));
    Workspace2D_sptr output2D =
        std::dynamic_pointer_cast<Workspace2D>(output_ws);

    // Test correct spectra loaded
    TS_ASSERT_EQUALS(output2D->getNumberHistograms(), specMax);
    // Check that spectra maps to correct detector
    auto detectorgroup = output2D->getSpectrum(specMax - 1).getDetectorIDs();
    TS_ASSERT_EQUALS(detectorgroup.size(), 1);
    TS_ASSERT_EQUALS(*detectorgroup.begin(), specMax);
  }

  void testLoadFailsIfEntryNumberOutOfRange() {
    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");
    ld.setPropertyValue("EntryNumber", "10");

    TS_ASSERT(!ld.isExecuted())
  }

  void testLoadFailsIfInvalidSpectraProperties() {
    std::vector<int> spectraIntegerList = {1, 123, 157};
    auto spectraList = createSpectraList(spectraIntegerList);
    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");
    ld.setPropertyValue("SpectrumList", spectraList);

    TS_ASSERT(!ld.isExecuted())
  }

  void testMaxThreadsRestoredWhenAlgorithmFinished() {
    int maxThreads = PARALLEL_GET_MAX_THREADS;
    LoadMuonNexusV2 ld;
    ld.initialize();
    ld.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    ld.setPropertyValue("OutputWorkspace", "outWS");

    ld.execute();

    TS_ASSERT_EQUALS(maxThreads, PARALLEL_GET_MAX_THREADS)
  }
};

//------------------------------------------------------------------------------
// Performance test
//------------------------------------------------------------------------------

class LoadMuonNexusV2TestPerformance : public CxxTest::TestSuite {
public:
  void setUp() override {
    loader.initialize();
    loader.setPropertyValue("Filename", "EMU00102347.nxs_v2");
    loader.setPropertyValue("OutputWorkspace", "ws");
  }

  void tearDown() override { AnalysisDataService::Instance().remove("ws"); }

  void testDefaultLoad() { loader.execute(); }

private:
  LoadMuonNexusV2 loader;
};