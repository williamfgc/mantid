#ifndef WORKSPACETEST_H_
#define WORKSPACETEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/SpectraDetectorMap.h"

using namespace Mantid::Kernel;
using namespace Mantid::API;

namespace Mantid { namespace DataObjects {
class WorkspaceTester : public MatrixWorkspace
{
public:
  WorkspaceTester() : MatrixWorkspace() {}
  virtual ~WorkspaceTester() {}

  // Empty overrides of virtual methods
  virtual const int getNumberHistograms() const { return 1;}
  const std::string id() const {return "WorkspaceTester";}
  void init(const int&, const int&, const int&)
  {
    // Put an 'empty' axis in to test the getAxis method
    m_axes.resize(1);
    m_axes[0] = new Axis(AxisType::Numeric,1);
  }
  int size() const {return 1;}
  int blocksize() const {return 1;}
  std::vector<double>& dataX(int const index) {return vec;}
  std::vector<double>& dataY(int const index) {return vec;}
  std::vector<double>& dataE(int const index) {return vec;}
  const std::vector<double>& dataX(int const index) const {return vec;}
  const std::vector<double>& dataY(int const index) const {return vec;}
  const std::vector<double>& dataE(int const index) const {return vec;}
  const IErrorHelper* errorHelper(int const index) const {return NULL;}
  void setErrorHelper(int const,IErrorHelper*) {}
  void setErrorHelper(int const,const IErrorHelper*) {}

private:
  std::vector<double> vec;
  int spec;
};
}} // namespace

DECLARE_WORKSPACE(WorkspaceTester)

class MatrixWorkspaceTest : public CxxTest::TestSuite
{
public:
  MatrixWorkspaceTest() : ws(new Mantid::DataObjects::WorkspaceTester)
  {
    ws->initialize(1,1,1);
  }
  
  void testGetSetTitle()
  {
    TS_ASSERT_EQUALS( ws->getTitle(), "" )
    ws->setTitle("something");
    TS_ASSERT_EQUALS( ws->getTitle(), "something" )
    ws->setTitle("");
  }

  void testGetSetComment()
  {
    TS_ASSERT_EQUALS( ws->getComment(), "" )
    ws->setComment("commenting");
    TS_ASSERT_EQUALS( ws->getComment(), "commenting" )
    ws->setComment("");
  }

  void testGetInstrument()
  {
    boost::shared_ptr<IInstrument> i = ws->getInstrument();
    TS_ASSERT_EQUALS( ws->getInstrument()->type(), "Instrument" )
  }

  void testSpectraMap()
  {
    MatrixWorkspace_sptr ws2 = WorkspaceFactory::Instance().create(ws,1,1,1);
    const SpectraDetectorMap &specs = ws2->spectraMap();
    TS_ASSERT_EQUALS( &(ws->spectraMap()), &specs )
    SpectraDetectorMap &specs2 = ws2->mutableSpectraMap();
    TS_ASSERT_DIFFERS( &(ws->spectraMap()), &specs2 )
  }

  void testGetSetSample()
  {
    TS_ASSERT( ws->getSample() )
    boost::shared_ptr<Sample> s(new Sample);
    TS_ASSERT_THROWS_NOTHING( ws->setSample(s) )
    TS_ASSERT_EQUALS( ws->getSample(), s )
    ws->getSample()->setName("test");
    TS_ASSERT_EQUALS( ws->getSample()->getName(), "test" )
  }

  void testGetMemorySize()
  {
    TS_ASSERT_THROWS_NOTHING( ws->getMemorySize() )
  }

  void testHistory()
  {
    TS_ASSERT_THROWS_NOTHING( WorkspaceHistory& h = ws->history() )
    const Mantid::DataObjects::WorkspaceTester wsc;
    TS_ASSERT_THROWS_NOTHING( const WorkspaceHistory& hh = wsc.getHistory() )
  }

  void testGetAxis()
  {
    TS_ASSERT_THROWS( ws->getAxis(-1), Exception::IndexError )
    TS_ASSERT_THROWS_NOTHING( ws->getAxis(0) )
    TS_ASSERT( ws->getAxis(0) )
    TS_ASSERT_THROWS( ws->getAxis(1), Exception::IndexError )
  }

  void testIsDistribution()
  {
    TS_ASSERT( ! ws->isDistribution() )
    TS_ASSERT( ws->isDistribution(true) )
    TS_ASSERT( ws->isDistribution() )
  }

  void testGetSetYUnit()
  {
    TS_ASSERT_EQUALS( ws->YUnit(), "Counts" )
    TS_ASSERT_THROWS_NOTHING( ws->setYUnit("something") )
    TS_ASSERT_EQUALS( ws->YUnit(), "something" )
  }

private:
  boost::shared_ptr<MatrixWorkspace> ws;

};

#endif /*WORKSPACETEST_H_*/
