// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_DATAHANDLING_DEFAULTEVENTLOADER_H_
#define MANTID_DATAHANDLING_DEFAULTEVENTLOADER_H_

#include "MantidAPI/Axis.h"
#include "MantidDataHandling/DllConfig.h"
#include "MantidDataHandling/EventWorkspaceCollection.h"

class BankPulseTimes;

namespace Mantid {
namespace DataHandling {
class LoadEventNexus;

/** Helper class for LoadEventNexus that is specific to the current default
  loading code for NXevent_data entries in Nexus files, in particular
  LoadBankFromDiskTask and ProcessBankData.
*/
class MANTID_DATAHANDLING_DLL DefaultEventLoader {
public:
  static void
  load(LoadEventNexus *alg, EventWorkspaceCollection &ws, bool haveWeights,
       bool event_id_is_spec, std::vector<std::string> bankNames,
       const std::vector<int> &periodLog, const std::string &classType,
       std::vector<std::size_t> bankNumEvents, const bool oldNeXusFileNames,
       const bool precount, const int chunk, const int totalChunks);

  /// Flag for dealing with a simulated file
  bool m_haveWeights;

  /// True if the event_id is spectrum no not pixel ID
  bool event_id_is_spec;

  /// whether or not to launch multiple ProcessBankData jobs per bank
  bool splitProcessing;

  /// Do we pre-count the # of events in each pixel ID?
  bool precount;

  /// Offset in the pixelID_to_wi_vector to use.
  detid_t pixelID_to_wi_offset;

  /// Maximum (inclusive) event ID possible for this instrument
  int32_t eventid_max{0};

  /// chunk number
  int chunk;
  /// number of chunks
  int totalChunks;
  /// for multiple chunks per bank
  int firstChunkForBank;
  /// number of chunks per bank
  size_t eventsPerChunk;

  LoadEventNexus *alg;
  EventWorkspaceCollection &m_ws;

  /// Vector where index = event_id; value = ptr to std::vector<TofEvent> in the
  /// event list.
  std::vector<std::vector<std::vector<Mantid::Types::Event::TofEvent> *>>
      eventVectors;

  /// linearize look up
  std::unordered_map<std::size_t, std::vector<Mantid::Types::Event::TofEvent> *>
      m_events;

  /// Vector where index = event_id; value = ptr to std::vector<WeightedEvent>
  /// in the event list.
  std::vector<std::vector<std::vector<Mantid::DataObjects::WeightedEvent> *>>
      weightedEventVectors;

  /// linearize look up
  std::unordered_map<std::size_t,
                     std::vector<Mantid::DataObjects::WeightedEvent> *>
      m_weightedEvents;

  /// Vector where (index = pixel ID+pixelID_to_wi_offset), value = workspace
  /// index)
  std::vector<size_t> pixelID_to_wi_vector;

  /// One entry of pulse times for each preprocessor
  std::vector<boost::shared_ptr<BankPulseTimes>> m_bankPulseTimes;

private:
  DefaultEventLoader(LoadEventNexus *alg, EventWorkspaceCollection &ws,
                     bool haveWeights, bool event_id_is_spec,
                     const size_t numBanks, const bool precount,
                     const int chunk, const int totalChunks);
  std::pair<size_t, size_t>
  setupChunking(std::vector<std::string> &bankNames,
                std::vector<std::size_t> &bankNumEvents);

  /// Sets event_max_id
  void SetEventIDMax();

  /// Map detector IDs to event lists.
  template <class T>
  void makeMapToEventLists(std::vector<std::vector<T>> &vectors);

  template <class T>
  void makeMapToEventLists(std::unordered_map<std::size_t, T> &umap);
};

template <class T>
void DefaultEventLoader::makeMapToEventLists(
    std::unordered_map<std::size_t, T> &umap) {

  // reserve maximum possible
  umap.reserve(m_ws.nPeriods() * (eventid_max + 1));

  if (event_id_is_spec) {

    for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
      for (size_t i = 0; i < m_ws.getNumberHistograms(); ++i) {
        const auto &spec = m_ws.getSpectrum(i);
        const size_t umapIndex =
            period * (eventid_max + 1) + spec.getSpectrumNo();
        getEventsFrom(m_ws.getSpectrum(i, period), umap[umapIndex]);
      }
    }

  } else {

    for (size_t j = 0; j < pixelID_to_wi_vector.size(); j++) {
      size_t wi = pixelID_to_wi_vector[j];
      // Save a POINTER to the vector
      if (wi < m_ws.getNumberHistograms()) {
        for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
          const size_t umapIndex =
              period * (eventid_max + 1) + (j - pixelID_to_wi_offset);
          getEventsFrom(m_ws.getSpectrum(wi, period), umap[umapIndex]);
        }
      }
    }
  }
}

/** Generate a look-up table where the index = the pixel ID of an event
 * and the value = a pointer to the EventList in the workspace
 * @param vectors :: the array to create the map on
 */
template <class T>
void DefaultEventLoader::makeMapToEventLists(
    std::vector<std::vector<T>> &vectors) {
  vectors.resize(m_ws.nPeriods());

  if (event_id_is_spec) {

    for (size_t i = 0; i < vectors.size(); ++i) {
      vectors[i].resize(eventid_max + 1, nullptr);
    }
    for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
      for (size_t i = 0; i < m_ws.getNumberHistograms(); ++i) {
        const auto &spec = m_ws.getSpectrum(i);
        getEventsFrom(m_ws.getSpectrum(i, period),
                      vectors[period][spec.getSpectrumNo()]);
      }
    }

  } else {
    // Make an array where index = pixel ID
    // Set the value to NULL by default
    for (size_t i = 0; i < vectors.size(); ++i) {
      vectors[i].resize(eventid_max + 1, nullptr);
    }

    for (size_t j = 0; j < pixelID_to_wi_vector.size(); j++) {
      size_t wi = pixelID_to_wi_vector[j];
      // Save a POINTER to the vector
      if (wi < m_ws.getNumberHistograms()) {
        for (size_t period = 0; period < m_ws.nPeriods(); ++period) {
          getEventsFrom(m_ws.getSpectrum(wi, period),
                        vectors[period][j - pixelID_to_wi_offset]);
        }
      }
    }
  }
}

} // namespace DataHandling
} // namespace Mantid

#endif /* MANTID_DATAHANDLING_DEFAULTEVENTLOADER_H_ */
