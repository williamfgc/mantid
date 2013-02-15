#include "MantidMDEvents/MDBox.h"
#include "MantidMDEvents/MDLeanEvent.h"
#include "MantidKernel/DiskBuffer.h"
#include "MantidMDEvents/MDGridBox.h"
#include "MantidMDEvents/BoxCtrlChangesList.h"

using Mantid::Kernel::DiskBuffer;
using namespace Mantid::API;

namespace Mantid
{
namespace MDEvents
{

//-----------------------------------------------------------------------------------------------
  /** Empty constructor */
  TMDE(MDBox)::MDBox()
   : MDBoxBase<MDE, nd>(),
     m_isLoaded(false),
     m_bIsMasked(false)
  {
  }

  //-----------------------------------------------------------------------------------------------
  /** Constructor
   * @param splitter :: BoxController that controls how boxes split
   * @param depth :: splitting depth of the new box.
   */
  TMDE(MDBox)::MDBox(BoxController_sptr splitter, const size_t depth,int64_t boxSize,int64_t boxID)
    : MDBoxBase<MDE, nd>(),
      m_isLoaded(false),
      m_bIsMasked(false)
  {
    if (splitter->getNDims() != nd)
      throw std::invalid_argument("MDBox::ctor(): controller passed has the wrong number of dimensions.");
    this->m_BoxController = splitter;
    this->m_depth = depth;

    if(boxID<0)  // Give it a fresh ID from the controller.
      this->setId( splitter->getNextId() );
    else         // somebody gives the ID on constructor
      this->setId(boxID);

    if(boxSize>0) data.reserve(boxSize);
   }

  //-----------------------------------------------------------------------------------------------
  /** Constructor
   * @param splitter :: BoxController that controls how boxes split
   * @param depth :: splitting depth of the new box.
   * @param extentsVector :: vector defining the extents
   */
  TMDE(MDBox)::MDBox(BoxController_sptr splitter, const size_t depth, const std::vector<Mantid::Geometry::MDDimensionExtents<coord_t> > & extentsVector,int64_t boxSize,int64_t boxID)
   :   MDBoxBase<MDE, nd>(extentsVector),
       m_isLoaded(false),
       m_bIsMasked(false)
  {
    if (splitter->getNDims() != nd)
      throw std::invalid_argument("MDBox::ctor(): controller passed has the wrong number of dimensions.");
    this->m_BoxController = splitter;
    this->m_depth = depth;
    // Give it a fresh ID from the controller.
    if(boxID<0) // Give it a fresh ID from the controller.
      this->setId( splitter->getNextId() );
    else     // somebody gives the ID on constructor
      this->setId(boxID);

    if(boxSize>0) data.reserve(boxSize);
  }


  //-----------------------------------------------------------------------------------------------
  /** Copy constructor
   * @param other: MDBox object to copy from.
   * */
  TMDE(MDBox)::MDBox(const MDBox<MDE,nd> & other)
     : MDBoxBase<MDE, nd>(other),
     data(other.data),
     m_isLoaded(other.m_isLoaded),
     m_bIsMasked(other.m_bIsMasked)
  {
  }


  //-----------------------------------------------------------------------------------------------
  /** Clear any points contained. */
  TMDE(
  void MDBox)::clear()
  {
    // Make sure the object is not in any of the disk MRUs, and mark any space it used as free
    //if (this->m_BoxController->useWriteBuffer())
    if(this->m_BoxController->isFileBacked())
       this->m_BoxController->getDiskBuffer().objectDeleted(this);
    // Clear all contents
    this->m_signal = 0.0;
    this->m_errorSquared = 0.0;

    this->clearDataFromMemory();

  }

  //-----------------------------------------------------------------------------------------------
  /** Clear the data[] vector ONLY but does not change the file-backed settings.
   * Used to free up the memory in a file-backed workspace without removing the events from disk. */
  TMDE(
  void MDBox)::clearDataFromMemory()
  {
    data.clear();
    vec_t().swap(data); // Linux trick to really free the memory
    m_isLoaded=false; 
    this->setBusy(false);
    // mark data unchanged
    this->resetDataChanges();

  }

  //-----------------------------------------------------------------------------------------------
  /** Returns the number of dimensions in this box */
  TMDE(
  size_t MDBox)::getNumDims() const
  {
    return nd;
  }

  //-----------------------------------------------------------------------------------------------
  /** Returns the number of un-split MDBoxes in this box (including all children)
   * @return :: 1 always since this is just a MD Box*/
  TMDE(
  size_t MDBox)::getNumMDBoxes() const
  {
    return 1;
  }

  //-----------------------------------------------------------------------------------------------
  /// Fill a vector with all the boxes up to a certain depth
  TMDE(
  void MDBox)::getBoxes(std::vector<MDBoxBase<MDE,nd> *> & boxes, size_t /*maxDepth*/, bool /*leafOnly*/)
  {
    boxes.push_back(this);
  }
  TMDE(
  void MDBox)::getBoxes(std::vector<Kernel::ISaveable *> & boxes, size_t /*maxDepth*/, bool /*leafOnly*/)
  {
    boxes.push_back(this);
  }


  //-----------------------------------------------------------------------------------------------
  /// Fill a vector with all the boxes up to a certain depth  
  TMDE(
  void MDBox)::getBoxes(std::vector<MDBoxBase<MDE,nd> *> & boxes, size_t /*maxDepth*/, bool /*leafOnly*/, Mantid::Geometry::MDImplicitFunction * /*function*/)
  {
    boxes.push_back(this);
  }
  TMDE(
  void MDBox)::getBoxes(std::vector<Kernel::ISaveable *> & boxes, size_t /*maxDepth*/, bool /*leafOnly*/, Mantid::Geometry::MDImplicitFunction * /*function*/)
  {
    boxes.push_back(this);
  }


  //-----------------------------------------------------------------------------------------------
  /** Returns the total number of points (events) in this box */
  TMDE(
  uint64_t MDBox)::getNPoints() const
  {
    if (this->wasSaved())
    {
      if (m_isLoaded)
        return data.size();
      else // m_fileNumEvents
        return this->getFileSize() + data.size();
    }
    else
      return data.size();
  }


 
  //-----------------------------------------------------------------------------------------------
  /** Returns a reference to the events vector contained within.
   * VERY IMPORTANT: call MDBox::releaseEvents() when you are done accessing that data.
   */
  TMDE(
  std::vector< MDE > & MDBox)::getEvents()
  {
    if (this->wasSaved())
    {
      // Load and concatenate the events if needed
      this->load();
      // The data vector is busy - can't release the memory yet
      this->setBusy();
      // Tell the to-write buffer to write out/discard the object (when no longer busy)
      this->m_BoxController->getDiskBuffer().toWrite(this);
    }
    // else: do nothing if the events are already in memory.
  // the non-const access to events assumes that the data will be modified;
    this->setDataChanged();
    return data;
  }

  TMDE(
  const std::vector<MDE> & MDBox)::getEvents()const 
  {
      return getConstEvents();
  }
  //-----------------------------------------------------------------------------------------------
  /** Returns a const reference to the events vector contained within.
   * VERY IMPORTANT: call MDBox::releaseEvents() when you are done accessing that data.
   */
  TMDE(
  const std::vector<MDE> & MDBox)::getConstEvents()const 
  {
    if (this->wasSaved())
    {
      // Load and concatenate the events if needed
      //TODO: redesighn
       MDBox<MDE,nd> *loader = const_cast<MDBox<MDE,nd> *>(this);
       loader->load();  // this will set isLoaded to true if not already loaded;
      // The data vector is busy - can't release the memory yet
      this->setBusy();

      // This access to data was const. Don't change the m_dataModified flag.

      // Tell the to-write buffer to discard the object (when no longer busy) as it has not been modified
      this->m_BoxController->getDiskBuffer().toWrite(this);
    }
    // else: do nothing if the events are already in memory.
    return data;
  }

  //-----------------------------------------------------------------------------------------------
  /** For file-backed MDBoxes, this marks that the data vector is
   * no longer "busy", and so it is safe for the MRU to cache it
   * back to disk if needed.
   */
  TMDE(
  void MDBox)::releaseEvents() 
  {
    // Data vector is no longer busy.
    this->setBusy(false);
  }


  //-----------------------------------------------------------------------------------------------
  /** Call to save the data (if needed) and release the memory used.
   *  Called from the DiskBuffer.
   *  If called directly presumes to know its file location and [TODO: refactor this] needs the file to be open correctly on correct group 
   */
  TMDE(
  void MDBox)::save()const
  {
  //      std::cout << "MDBox ID " << this->getId() << " being saved." << std::endl;


   // this aslo indirectly checks if the object knows its place (may be wrong place but no checks for that here)
   if (this->wasSaved())
   {
     //TODO: redesighn const_cast
    // This will load and append events ONLY if needed.
      MDBox<MDE,nd> *loader = const_cast<MDBox<MDE,nd> *>(this);
      loader->load();  // this will set isLoaded to true if not already loaded;
   

      // This is the new size of the event list, possibly appended (if used AddEvent) or changed otherwise (non-const access)
      if (data.size() > 0)
      {

         // Save at the ISaveable specified place
          this->saveNexus(this->m_BoxController->getFile());
      }
   } 
   else
     if(data.size()>0)  throw std::runtime_error(" Attempt to save undefined event");
   
   
  }


  //-----------------------------------------------------------------------------------------------
  /** Save the box's Event data to an open nexus file.
   *
   * @param file :: Nexus File object, must already by opened with MDE::prepareNexusData()
   */
  TMDE(
  inline void MDBox)::saveNexus(::NeXus::File * file) const
  {
    //std::cout << "Box " << this->getId() << " saving to " << m_fileIndexStart << std::endl;
    MDE::saveVectorToNexusSlab(this->data, file, this->getFilePosition(),
                               this->m_signal, this->m_errorSquared);

  }


  //-----------------------------------------------------------------------------------------------
  /** Load the box's Event data from an open nexus file.
   * The FileIndex start and numEvents must be set correctly already.
   * Clear existing data from memory!
   *
   * @param file :: Nexus File object, must already by opened with MDE::openNexusData()
   */
  TMDE(
  inline void MDBox)::loadNexus(::NeXus::File * file, bool setIsLoaded)
  {
    this->data.clear();
    uint64_t fileIndexStart = this->getFilePosition();
    uint64_t fileNumEvents  = this->getFileSize();
    if(fileIndexStart == std::numeric_limits<uint64_t>::max())
      throw(std::runtime_error("MDBox::loadNexus -- attempt to load box from undefined location"));
    MDE::loadVectorFromNexusSlab(this->data, file, fileIndexStart, fileNumEvents);

   
    this->m_isLoaded=setIsLoaded;
  
  }


 TMDE(
 inline void MDBox)::load()
 {
    // Is the data in memory right now (cached copy)?
    if (!m_isLoaded)
    {
      // Perform the data loading
      ::NeXus::File * file = this->m_BoxController->getFile();
      if (file)
      {
        // Mutex for disk access (prevent read/write at the same time)
        Kernel::RecursiveMutex & mutex = this->m_BoxController->getDiskBuffer().getFileMutex();
        mutex.lock();
        // Note that this APPENDS any events to the existing event list
        //  (in the event that addEvent() was called for a box that was on disk)
        try
        {
          uint64_t fileIndexStart = this->getFilePosition();
          uint64_t fileNumEvents  = this->getFileSize();
          MDE::loadVectorFromNexusSlab(data, file,fileIndexStart, fileNumEvents);
          m_isLoaded = true;
          mutex.unlock();
        }
        catch (std::exception &e)
        {
          mutex.unlock();
          throw e;
        }
      }
    }
  }



  //-----------------------------------------------------------------------------------------------
  /** Allocate and return a vector with a copy of all events contained
   */
  TMDE(
  std::vector< MDE > * MDBox)::getEventsCopy()
  {
    std::vector< MDE > * out = new std::vector<MDE>();
    //Make the copy
    out->insert(out->begin(), data.begin(), data.end());
    return out;
  }

  //-----------------------------------------------------------------------------------------------
  /** Refresh the cache.
   *
   * For MDBox, if MDBOX_TRACK_SIGNAL_WHEN_ADDING is defined,
   * then the signal/error is tracked on adding, so
   * this does nothing.

   * beware, that it wrongly accumulates signal and error when part of the data is on file and 
   * and some recent data were not saved to the HDD before adding new data
   * This actually means that refreshCache has to be called only once in events adding process
   */
  TMDE(
  void MDBox)::refreshCache(Kernel::ThreadScheduler * /*ts*/)
  {
#ifndef MDBOX_TRACK_SIGNAL_WHEN_ADDING
    // Use the cached value if it is on disk
    double signalSum(0);
    double errorSum(0);

    if (this->wasSaved()) // There are possible problems with disk buffered events, as saving calculates averages and these averages has to be added to memory contents
    {
      if(!m_isLoaded)  // events were saved,  averages calculated and stored 
      {
        // the partial data were not loaded from HDD but their averages should be calculated when loaded. Add them 
         signalSum +=double(this->m_signal);
         errorSum  +=double(this->m_errorSquared);
      }
    }
    // calculate all averages from memory
    typename std::vector<MDE>::const_iterator it_end = data.end();
    for(typename std::vector<MDE>::const_iterator it = data.begin(); it != it_end; ++it)
    {
        const MDE & event = *it;
        // Convert floats to doubles to preserve precision when adding them.
        signalSum += static_cast<signal_t>(event.getSignal());
        errorSum += static_cast<signal_t>(event.getErrorSquared());
    }

    this->m_signal = signal_t(signalSum);
    this->m_errorSquared=signal_t(errorSum);

    /// TODO #4734: sum the individual weights of each event?
    this->m_totalWeight = static_cast<double>(this->getNPoints());
#endif
  }


  //-----------------------------------------------------------------------------------------------
  /** Calculate and ccache the centroid of this box.
   */
  TMDE(
  void MDBox)::refreshCentroid(Kernel::ThreadScheduler * /*ts*/)
  {
#ifdef MDBOX_TRACK_CENTROID
    calculateCentroid(this->m_centroid);
#endif
  }

  //-----------------------------------------------------------------------------------------------
  /** Calculate the centroid of this box.
   * @param centroid [out] :: nd-sized array that will be set to the centroid.
   */
  TMDE(
  void MDBox)::calculateCentroid(coord_t * centroid) const
  {
    for (size_t d=0; d<nd; d++)
      centroid[d] = 0;

    // Signal was calculated before (when adding)
    // Keep 0.0 if the signal is null. This avoids dividing by 0.0
    if (this->m_signal == 0) return;

    typename std::vector<MDE>::const_iterator it_end = data.end();
    for(typename std::vector<MDE>::const_iterator it = data.begin(); it != it_end; ++it)
    {
      const MDE & Evnt = *it;
      double signal = Evnt.getSignal();
      for (size_t d=0; d<nd; d++)
      {
        // Total up the coordinate weighted by the signal.
        centroid[d] += Evnt.getCenter(d) * static_cast<coord_t>(signal);
      }
    }

    // Normalize by the total signal
    for (size_t d=0; d<nd; d++)
      centroid[d] /= coord_t(this->m_signal);
  }

  //-----------------------------------------------------------------------------------------------
  /** Calculate the statistics for each dimension of this MDBox, using
   * all the contained events
   * @param stats :: nd-sized fixed array of MDDimensionStats, reset to 0.0 before!
   */
  TMDE(
  void MDBox)::calculateDimensionStats(MDDimensionStats * stats) const
  {
    typename std::vector<MDE>::const_iterator it_end = data.end();
    for(typename std::vector<MDE>::const_iterator it = data.begin(); it != it_end; ++it)
    {
      const MDE & Evnt = *it;
      for (size_t d=0; d<nd; d++)
      {
        stats[d].addPoint( Evnt.getCenter(d) );
      }
    }
  }

  //-----------------------------------------------------------------------------------------------
  /** Add a MDLeanEvent to the box.
   * @param event :: reference to a MDEvent to add.
   * */
  TMDE(
  void MDBox)::addEvent( const MDE & Evnt)
  {
    dataMutex.lock();
    this->data.push_back(Evnt );

//    // When we reach the split threshold exactly, track that the MDBox is too small
//    // We check on equality and not >= to only add a box once.
//    if (this->data.size() == this->m_BoxController->getSplitThreshold())
//    {
//      this->m_BoxController->addBoxToSplit(this);
//    }

#ifdef MDBOX_TRACK_SIGNAL_WHEN_ADDING
    // Keep the running total of signal and error
    signal_t signal = static_cast<signal_t>(Evnt.getSignal());
    this->m_signal += signal;
    this->m_errorSquared += static_cast<signal_t>(Evnt.getErrorSquared());
#endif

#ifdef MDBOX_TRACKCENTROID_WHENADDING
    // Running total of the centroid
    for (size_t d=0; d<nd; d++)
      this->m_centroid[d] += Evnt.getCenter(d) * signal;
#endif

    dataMutex.unlock();
  }

  //-----------------------------------------------------------------------------------------------
  /** Add a MDEvent to the box.
   // add a single event and set pounter to the box which needs splitting (if one actually need) 

   * @param point :: reference to a MDEvent to add.
   * @returns  :: pointer to itself if the box should be split or NULL if not yet
   * */
  TMDE(
  void MDBox)::addAndTraceEvent( const MDE & point, size_t index)
  {
    dataMutex.lock();
    this->data.push_back(point);

    dataMutex.unlock();

    // When we reach the split threshold exactly, track that the MDBox is too small
    // We check on equality and not >= to only add a box once.
    if (this->data.size() == this->m_BoxController->getSplitThreshold())
    {     
       auto BoxCtrl = dynamic_cast<BoxCtrlChangesList<MDBoxToChange<MDE,nd> >*>(this->m_BoxController.get());
       BoxCtrl->addBoxToSplit(MDBoxToChange<MDE,nd>(this,index));

    }



  }

  //-----------------------------------------------------------------------------------------------
  /** Add a MDLeanEvent to the box, in a NON-THREAD-SAFE manner.
   * No lock is performed. This is only safe if no 2 threads will
   * try to add to the same box at the same time.
   *
   * @param event :: reference to a MDEvent to add.
   * */
  TMDE(
  void MDBox)::addEventUnsafe( const MDE & Evnt)
  {
    this->data.push_back(Evnt );


#ifdef MDBOX_TRACK_SIGNAL_WHEN_ADDING
    // Keep the running total of signal and error
    double signal = static_cast<signal_t>(Evnt.getSignal());
    this->m_signal += signal;
    this->m_errorSquared += static_cast<signal_t>(Evnt.getErrorSquared());
#endif

#ifdef MDBOX_TRACKCENTROID_WHENADDING
    // Running total of the centroid
    for (size_t d=0; d<nd; d++)
      this->m_centroid[d] += Evnt.getCenter(d) * signal;
#endif
  }

  //-----------------------------------------------------------------------------------------------
  /** Add several events. No bounds checking is made!
   *
   * @param events :: vector of events to be copied.
   * @param start_at :: index to start at in vector
   * @param stop_at :: index to stop at in vector (exclusive)
   * @return the number of events that were rejected (because of being out of bounds)
   */
  TMDE(
  size_t MDBox)::addEventsPart(const std::vector<MDE> & events, const size_t start_at, const size_t stop_at)
  {
    dataMutex.lock();
    typename std::vector<MDE>::const_iterator start = events.begin()+start_at;
    typename std::vector<MDE>::const_iterator end = events.begin()+stop_at;
    // Copy all the events
    this->data.insert(this->data.end(), start, end);

 
#ifdef MDBOX_TRACK_SIGNAL_WHEN_ADDING
    //Running total of signal/error
    for(typename std::vector<MDE>::const_iterator it = start; it != end; ++it)
    {
      double signal = static_cast<signal_t>(it->getSignal());
      this->m_signal += signal;
      this->m_errorSquared += static_cast<signal_t>(it->getErrorSquared());

#ifdef MDBOX_TRACKCENTROID_WHENADDING
    // Running total of the centroid
    for (size_t d=0; d<nd; d++)
      this->m_centroid[d] += it->getCenter(d) * signal;
#endif
    }
#endif

    dataMutex.unlock();
    return 0;
  }


  //-----------------------------------------------------------------------------------------------
  /** Add several events, within a given range, with no bounds checking,
   * and not in a thread-safe way
   *
   * @param events :: vector of events to be copied.
   * @param start_at :: index to start at in vector
   * @param stop_at :: index to stop at in vector (exclusive)
   * @return the number of events that were rejected (because of being out of bounds)
   */
  TMDE(
  size_t MDBox)::addEventsPartUnsafe(const std::vector<MDE> & events, const size_t start_at, const size_t stop_at)
  {
    // The regular MDBox is just as safe/unsafe
    return this->addEventsPart(events, start_at, stop_at);
  }


  //-----------------------------------------------------------------------------------------------
  /** Perform centerpoint binning of events.
   * @param bin :: MDBin object giving the limits of events to accept.
   * @param fullyContained :: optional bool array sized [nd] of which dimensions are known to be fully contained (for MDSplitBox)
   */
  TMDE(
  void MDBox)::centerpointBin(MDBin<MDE,nd> & bin, bool * fullyContained) const
  {
    if (fullyContained)
    {
      // For MDSplitBox, check if we've already found that all dimensions are fully contained
      size_t d;
      for (d=0; d<nd; ++d)
      {
        if (!fullyContained[d]) break;
      }
      if (d == nd)
      {
        // All dimensions are fully contained, so just return the cached total signal instead of counting.
        bin.m_signal += static_cast<signal_t>(this->m_signal);
        bin.m_errorSquared += static_cast<signal_t>(this->m_errorSquared);
        return;
      }
    }

    // If the box is cached to disk, you need to retrieve it
    const std::vector<MDE> & events = this->getConstEvents();
    typename std::vector<MDE>::const_iterator it = events.begin();
    typename std::vector<MDE>::const_iterator it_end = events.end();

    // For each MDLeanEvent
    for (; it != it_end; ++it)
    {
      // Go through each dimension
      size_t d;
      for (d=0; d<nd; ++d)
      {
        // Check that the value is within the bounds given. (Rotation is for later)
        coord_t x = it->getCenter(d);
        if (x < bin.m_min[d])
          break;
        if (x >= bin.m_max[d])
          break;
      }
      // If the loop reached the end, then it was all within bounds.
      if (d == nd)
      {
        // Accumulate error and signal (as doubles, to preserve precision)
        bin.m_signal += static_cast<signal_t>(it->getSignal());
        bin.m_errorSquared += static_cast<signal_t>(it->getErrorSquared());
      }
    }
    // it is constant access, so no saving or fiddling with the buffer is needed. Events just can be dropped if necessary
    this->setBusy(false);
    //this->releaseEvents();
  }


  //-----------------------------------------------------------------------------------------------
  /** General (non-axis-aligned) centerpoint binning method.
   * TODO: TEST THIS!
   *
   * @param bin :: a MDBin object giving the limits, aligned with the axes of the workspace,
   *        of where the non-aligned bin MIGHT be present.
   * @param function :: a ImplicitFunction that will evaluate true for any coordinate that is
   *        contained within the (non-axis-aligned) bin.
   */
  TMDE(
  void MDBox)::generalBin(MDBin<MDE,nd> & bin, Mantid::Geometry::MDImplicitFunction & function) const
  {
    UNUSED_ARG(bin);

    typename std::vector<MDE>::const_iterator it = data.begin();
    typename std::vector<MDE>::const_iterator it_end = data.end();
    // For each MDLeanEvent
    for (; it != it_end; ++it)
    {
      if (function.isPointContained(it->getCenter())) //HACK
      {
        // Accumulate error and signal
        bin.m_signal += static_cast<signal_t>(it->getSignal());
        bin.m_errorSquared += static_cast<signal_t>(it->getErrorSquared());
      }
    }
  }


  /** Integrate the signal within a sphere; for example, to perform single-crystal
   * peak integration.
   * The CoordTransform object could be used for more complex shapes, e.g. "lentil" integration, as long
   * as it reduces the dimensions to a single value.
   *
   * @param radiusTransform :: nd-to-1 coordinate transformation that converts from these
   *        dimensions to the distance (squared) from the center of the sphere.
   * @param radiusSquared :: radius^2 below which to integrate
   * @param[out] signal :: set to the integrated signal
   * @param[out] errorSquared :: set to the integrated squared error.
   */
  TMDE(
  void MDBox)::integrateSphere(Mantid::API::CoordTransform & radiusTransform, const coord_t radiusSquared, signal_t & signal, signal_t & errorSquared) const
  {
    // If the box is cached to disk, you need to retrieve it
    const std::vector<MDE> & events = this->getConstEvents();
    typename std::vector<MDE>::const_iterator it = events.begin();
    typename std::vector<MDE>::const_iterator it_end = events.end();

    // For each MDLeanEvent
    for (; it != it_end; ++it)
    {
      coord_t out[nd];
      radiusTransform.apply(it->getCenter(), out);
      if (out[0] < radiusSquared)
      {
        signal += static_cast<signal_t>(it->getSignal());
        errorSquared += static_cast<signal_t>(it->getErrorSquared());
      }
    }
    // it is constant access, so no saving or fiddling with the buffer is needed. Events just can be dropped if necessary
    this->setBusy(false);
    //this->releaseEvents();
  }

  //-----------------------------------------------------------------------------------------------
  /** Find the centroid of all events contained within by doing a weighted average
   * of their coordinates.
   *
   * @param radiusTransform :: nd-to-1 coordinate transformation that converts from these
   *        dimensions to the distance (squared) from the center of the sphere.
   * @param radiusSquared :: radius^2 below which to integrate
   * @param[out] centroid :: array of size [nd]; its centroid will be added
   * @param[out] signal :: set to the integrated signal
   */
  TMDE(
  void MDBox)::centroidSphere(Mantid::API::CoordTransform & radiusTransform, const coord_t radiusSquared, coord_t * centroid, signal_t & signal) const
  {
    // If the box is cached to disk, you need to retrieve it
    const std::vector<MDE> & events = this->getConstEvents();
    typename std::vector<MDE>::const_iterator it = events.begin();
    typename std::vector<MDE>::const_iterator it_end = events.end();

    // For each MDLeanEvent
    for (; it != it_end; ++it)
    {
      coord_t out[nd];
      radiusTransform.apply(it->getCenter(), out);
      if (out[0] < radiusSquared)
      {
        coord_t eventSignal = static_cast<coord_t>(it->getSignal());
        signal += signal_t(eventSignal);
        for (size_t d=0; d<nd; d++)
          centroid[d] += it->getCenter(d) * eventSignal;
      }
    }
    // it is constant access, so no saving or fiddling with the buffer is needed. Events just can be dropped if necessary
    this->setBusy(false);
  }


  //-----------------------------------------------------------------------------------------------
  /** Transform the dimensions contained in this box
   * x' = x*scaling + offset.
   *
   * @param scaling :: multiply each coordinate by this value.
   * @param offset :: after multiplying, add this offset.
   */
  TMDE(
  void MDBox)::transformDimensions(std::vector<double> & scaling, std::vector<double> & offset)
  {
    MDBoxBase<MDE,nd>::transformDimensions(scaling, offset);
    std::vector<MDE> & events = this->getEvents();
    typename std::vector<MDE>::iterator it;
    typename std::vector<MDE>::iterator it_end = events.end();
    for (it = events.begin(); it != it_end; ++it)
    {
      coord_t * center = it->getCenterNonConst();
      for (size_t d=0; d<nd; d++)
        center[d] = (center[d] * static_cast<coord_t>(scaling[d])) + static_cast<coord_t>(offset[d]);
    }
    this->setBusy(false);
  }

    ///Setter for masking the box
  TMDE(
  void MDBox)::mask()
  {
    m_bIsMasked = true;
  }

  ///Setter for unmasking the box
  TMDE(
  void MDBox)::unmask()
  {
    m_bIsMasked = false;
  }



}//namespace MDEvents

}//namespace Mantid

