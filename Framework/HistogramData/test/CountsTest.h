#ifndef MANTID_HISTOGRAMDATA_COUNTSTEST_H_
#define MANTID_HISTOGRAMDATA_COUNTSTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidHistogramData/BinEdges.h"
#include "MantidHistogramData/Counts.h"
#include "MantidHistogramData/Frequencies.h"

using Mantid::HistogramData::BinEdges;
using Mantid::HistogramData::Counts;
using Mantid::HistogramData::Frequencies;

class CountsTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static CountsTest *createSuite() { return new CountsTest(); }
  static void destroySuite(CountsTest *suite) { delete suite; }

  void test_construct_default() {
    const Counts counts{};
    TS_ASSERT(!counts);
  }

  void test_construct_from_null_Frequencies() {
    const Frequencies frequencies{};
    const BinEdges edges{};
    const Counts counts(frequencies, edges);
    TS_ASSERT(!counts);
  }

  void test_construct_from_empty_Frequencies() {
    const Frequencies frequencies(0);
    const BinEdges edges{0};
    const Counts counts(frequencies, edges);
    TS_ASSERT_EQUALS(counts.size(), 0);
  }

  void test_construct_from_empty_Frequencies_null_BinEdges() {
    const Frequencies frequencies(0);
    const BinEdges edges{};
    TS_ASSERT_THROWS(const Counts counts(frequencies, edges), std::logic_error);
  }

  void test_construct_from_empty_Frequencies_size_mismatch() {
    const Frequencies frequencies(0);
    const BinEdges edges(2);
    TS_ASSERT_THROWS(const Counts counts(frequencies, edges), std::logic_error);
  }

  void test_construct_from_Frequencies_null_BinEdges() {
    const Frequencies frequencies(1);
    const BinEdges edges{};
    TS_ASSERT_THROWS(const Counts counts(frequencies, edges), std::logic_error);
  }

  void test_construct_from_Frequencies_size_mismatch() {
    const Frequencies frequencies(2);
    const BinEdges edges(2);
    TS_ASSERT_THROWS(const Counts counts(frequencies, edges), std::logic_error);
  }

  void test_construct_from_Frequencies() {
    const Frequencies frequencies{1.0, 2.0};
    const BinEdges edges{0.1, 0.2, 0.4};
    const Counts counts(frequencies, edges);
    TS_ASSERT_EQUALS(counts.size(), 2);
    TS_ASSERT_DELTA(counts[0], 0.1, 1e-14);
    TS_ASSERT_DELTA(counts[1], 0.4, 1e-14);
  }

  void test_move_construct_from_Frequencies() {
    Frequencies frequencies(1);
    const BinEdges edges(2);
    auto old_ptr = &frequencies[0];
    const Counts counts(std::move(frequencies), edges);
    TS_ASSERT(!frequencies);
    TS_ASSERT_EQUALS(&counts[0], old_ptr);
  }

  void test_move_construct_from_Frequencies_and_cow() {
    Frequencies frequencies(1);
    const Frequencies copy(frequencies);
    const BinEdges edges(2);
    auto old_ptr = &frequencies[0];
    const Counts counts(std::move(frequencies), edges);
    // Moved from frequencies...
    TS_ASSERT(!frequencies);
    // ... but made a copy of data, since "copy" also held a reference.
    TS_ASSERT_DIFFERS(&counts[0], old_ptr);
  }
};

#endif /* MANTID_HISTOGRAMDATA_COUNTSTEST_H_ */
