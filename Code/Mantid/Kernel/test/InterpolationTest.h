#ifndef INTERPOLATIONTEST_H_
#define INTERPOLATIONTEST_H_

#include <cxxtest/TestSuite.h>
#include <ctime>
#include "MantidKernel/Interpolation.h"

using namespace Mantid::Kernel;

class InterpolationTest : public CxxTest::TestSuite
{
public:
  
	void test1()
	{
    Interpolation interpolation;

    interpolation.addPoint(200.0, 50);
    interpolation.addPoint(201.0, 60);
    interpolation.addPoint(202.0, 100);
    interpolation.addPoint(204.0, 400);
    interpolation.addPoint(203.0, 300);

    // Test that all the base class member variables are correctly assigned to
    TS_ASSERT_DELTA( interpolation.value(100), 50.0 ,0.000000001); 
    TS_ASSERT_DELTA( interpolation.value(3000), 400.0 ,0.000000001);
    TS_ASSERT_DELTA( interpolation.value(200.5), 55.0 ,0.000000001); 
    TS_ASSERT_DELTA( interpolation.value(201.25), 70.0 ,0.000000001); 
    TS_ASSERT_DELTA( interpolation.value(203.5), 350.0 ,0.000000001); 


    interpolation.setXUnit("bob");
    std::stringstream str;
    str << interpolation;
    TS_ASSERT( str.str().compare("linear ; bob ; 200 50 ; 201 60 ; 202 100 ; 203 300 ; 204 400") == 0 );

    Interpolation readIn;
    TS_ASSERT( readIn.getXUnit().compare("TOF") == 0 );
    str >> readIn;
    TS_ASSERT( readIn.getXUnit().compare("bob") == 0 );

    // Test that all the base class member variables are correctly assigned to
    TS_ASSERT_DELTA( readIn.value(100), 50.0 ,0.000000001); 
    TS_ASSERT_DELTA( readIn.value(3000), 400.0 ,0.000000001);
    TS_ASSERT_DELTA( readIn.value(200.5), 55.0 ,0.000000001); 
    TS_ASSERT_DELTA( readIn.value(201.25), 70.0 ,0.000000001); 
    TS_ASSERT_DELTA( readIn.value(203.5), 350.0 ,0.000000001); 
	}

	void testEmpty()
	{
    Interpolation interpolation;

    std::stringstream str;
    str << interpolation;
    TS_ASSERT( str.str().compare("linear ; TOF") == 0 );

    Interpolation readIn;
    str >> readIn;

    TS_ASSERT( readIn.containData() == false );
	}
  
};

#endif /*INTERPOLATIONTEST_H_*/
