#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/LocatedDataRef.h"
#include "MantidAPI/WorkspaceIterator.h"
#include "MantidAPI/WorkspaceIteratorCode.h"
#include "MantidAPI/SpectraDetectorMap.h"
#include "MantidAPI/ParInstrument.h"
#include "MantidGeometry/DetectorGroup.h"

namespace Mantid
{
namespace API
{

Kernel::Logger& MatrixWorkspace::g_log = Kernel::Logger::get("Workspace");

/// Default constructor
MatrixWorkspace::MatrixWorkspace() : Workspace(), m_axes(), m_isInitialized(false),
  sptr_instrument(new Instrument), m_spectramap(), sptr_sample(new Sample),
   m_YUnit("Counts"), m_isDistribution(false),sptr_parmap(new Geometry::ParameterMap)
{}

/// Destructor
// RJT, 3/10/07: The Analysis Data Service needs to be able to delete workspaces, so I moved this from protected to public.
MatrixWorkspace::~MatrixWorkspace()
{
  for (unsigned int i = 0; i < m_axes.size(); ++i)
  {
    delete m_axes[i];
  }
}

/** Initialize the workspace. Calls the protected init() method, which is implemented in each type of
 *  workspace. Returns immediately if the workspace is already initialized.
 *  @param NVectors The number of spectra in the workspace (only relevant for a 2D workspace
 *  @param XLength The number of X data points/bin boundaries in each vector (must all be the same)
 *  @param YLength The number of data/error points in each vector (must all be the same)
 */
void MatrixWorkspace::initialize(const int &NVectors, const int &XLength, const int &YLength)
{
  // Check validity of arguments
  if (NVectors <= 0 || XLength <= 0 || YLength <= 0)
  {
    g_log.error("All arguments to init must be positive and non-zero");
    throw std::out_of_range("All arguments to init must be positive and non-zero");
  }

  // Bypass the initialization if the workspace has already been initialized.
  if (m_isInitialized) return;

  // Invoke init() method of the derived class inside a try/catch clause
  try
  {
    this->init(NVectors, XLength, YLength);
  }
  catch(std::runtime_error& ex)
  {
    g_log.error() << "Error initializing the workspace" << ex.what() << std::endl;
    throw;
  }

  // Indicate that this Algorithm has been initialized to prevent duplicate attempts.
  m_isInitialized = true;
}

/** Set the instrument
 *
 * \param instr Shared pointer to an instrument.
 */
void MatrixWorkspace::setInstrument(const IInstrument_sptr& instr)
{
    boost::shared_ptr<Instrument> tmp = boost::dynamic_pointer_cast<Instrument>(instr);
    if (tmp)
    {
        sptr_instrument=tmp;
    }
    else
    {
        boost::shared_ptr<ParInstrument> tmp = boost::dynamic_pointer_cast<ParInstrument>(instr);
        if (tmp)
        {
            sptr_instrument = tmp->baseInstrument();
            sptr_parmap = tmp->getParameterMap();
        }
    }
        
}

/** Set the sample
 *
 *  @param sample A shared pointer to the sample
 */
void MatrixWorkspace::setSample(const boost::shared_ptr<Sample>& sample)
{
  sptr_sample = sample;
}

/** Get a const reference to the SpectraDetectorMap associated with this workspace.
 *  Can ONLY be taken as a const reference!
 *
 *  @return The SpectraDetectorMap
 */
const SpectraDetectorMap& MatrixWorkspace::spectraMap() const
{
  return *m_spectramap;
}

/** Get a reference to the SpectraDetectorMap associated with this workspace.
 *  This non-const method will copy the map if it is shared between more than one workspace,
 *  and the reference returned will be to the copy.
 *  Can ONLY be taken by reference!
 *
 *  @return The SpectraDetectorMap
 */
SpectraDetectorMap& MatrixWorkspace::mutableSpectraMap()
{
  return m_spectramap.access();
}

/** Get the effective detector for the given spectrum
 *  @param  index The workspace index for which the detector is required
 *  @return A single detector object representing the detector(s) contributing
 *          to the given spectrum number. If more than one detector contributes then
 *          the returned object's concrete type will be DetectorGroup.
 *  @throw  std::runtime_error if the SpectraDetectorMap has not been filled
 */
Geometry::IDetector_sptr MatrixWorkspace::getDetector(const int index) const
{
  if ( ! m_spectramap->nElements() )
  {
    g_log.error("SpectraDetectorMap has not been populated.");
    throw std::runtime_error("SpectraDetectorMap has not been populated.");
  }
  
  const int spectrum_number = getAxis(1)->spectraNo(index);
  const std::vector<int> dets = m_spectramap->getDetectors(spectrum_number);
  if ( dets.empty() )
  {
    g_log.debug() << "Spectrum number " << spectrum_number << " not found" << std::endl;
    throw Kernel::Exception::NotFoundError("Spectrum number not found", spectrum_number);
  }
  else if ( dets.size() == 1 ) 
  {
    // If only 1 detector for the spectrum number, just return it
    return getInstrument()->getDetector(dets[0]);
  }
  // Else need to construct a DetectorGroup and return that
  std::vector<Geometry::IDetector_sptr> dets_ptr;
  std::vector<int>::const_iterator it;
  for ( it = dets.begin(); it != dets.end(); ++it )
  {
    dets_ptr.push_back( getInstrument()->getDetector(*it) );
  }
  
  return Geometry::IDetector_sptr( new Geometry::DetectorGroup(dets_ptr) );
}

/** Returns the 2Theta scattering angle for a detector
 *  @param det A pointer to the detector object (N.B. might be a DetectorGroup)
 *  @return The scattering angle (0 < theta < pi)
 */
const double MatrixWorkspace::detectorTwoTheta(Geometry::IDetector_const_sptr det) const
{
  const Geometry::V3D samplePos = getInstrument()->getSample()->getPos();
  const Geometry::V3D beamLine = samplePos - getInstrument()->getSource()->getPos();
  const Geometry::V3D sampleDetVec = det->getPos() - samplePos;
  return sampleDetVec.angle(beamLine);
}

/** Get a shared pointer to the instrument associated with this workspace
 *
 *  @return The instrument class
 */
IInstrument_sptr MatrixWorkspace::getInstrument()const
{
    if ( sptr_parmap->size() == 0 )  return sptr_instrument;
    ParInstrument* pi = new ParInstrument(sptr_instrument,sptr_parmap);
    IInstrument* ii = static_cast<IInstrument*>(pi);
    boost::shared_ptr<IInstrument> tmp(ii);
    return (tmp);
}

/** Get a shared pointer to the instrument associated with this workspace
 *
 *  @return The instrument class
 */
boost::shared_ptr<Instrument> MatrixWorkspace::getBaseInstrument()const
{
    return sptr_instrument;
}

/** Get the sample associated with this workspace
 *
 *  @return The sample class
 */
boost::shared_ptr<Sample> MatrixWorkspace::getSample() const
{
  return sptr_sample;
}

/// The number of axes which this workspace has
const int MatrixWorkspace::axes() const
{
  return static_cast<int>(m_axes.size());
}

/** Get a pointer to a workspace axis
 *  @param axisIndex The index of the axis required
 *  @throw IndexError If the argument given is outside the range of axes held by this workspace
 */
Axis* const MatrixWorkspace::getAxis(const int axisIndex) const
{
  if ( axisIndex < 0 || axisIndex >= static_cast<int>(m_axes.size()) )
  {
    g_log.error() << "Argument to getAxis (" << axisIndex << ") is invalid for this (" << m_axes.size() << " axis) workspace" << std::endl;
    throw Kernel::Exception::IndexError(axisIndex, m_axes.size(),"Argument to getAxis is invalid for this workspace");
  }

  return m_axes[axisIndex];
}

/// Returns true if the workspace contains data in histogram form, false if it's point-like
const bool MatrixWorkspace::isHistogramData() const
{
  return ( readX(0).size()==readY(0).size() ? false : true );
}

/// Returns the units of the data in the workspace (default 'Counts')
std::string MatrixWorkspace::YUnit() const
{
  std::string retVal = m_YUnit;
  // If this workspace a distribution & has at least one axis & this axis has its unit set
  // then append that unit to the string to be returned
  if ( this->isDistribution() && this->axes() && this->getAxis(0)->unit() )
  {
    retVal = retVal + " per " + this->getAxis(0)->unit()->label();
  }
  return retVal;
}

/// Sets a new unit for the data (Y axis) in the workspace
void MatrixWorkspace::setYUnit(const std::string& newUnit)
{
  m_YUnit = newUnit;
}

/// Are the Y-values in this workspace dimensioned?
const bool& MatrixWorkspace::isDistribution() const
{
  return m_isDistribution;
}

/// Set the flag for whether the Y-values are dimensioned
bool& MatrixWorkspace::isDistribution(bool newValue)
{
  m_isDistribution = newValue;
  return m_isDistribution;
}

long int MatrixWorkspace::getMemorySize() const
{
    return 3*size()*sizeof(double)/1024;
}


} // namespace API
} // Namespace Mantid


///\cond TEMPLATE
template DLLExport class Mantid::API::workspace_iterator<Mantid::API::LocatedDataRef,Mantid::API::MatrixWorkspace>;
template DLLExport class Mantid::API::workspace_iterator<const Mantid::API::LocatedDataRef, const Mantid::API::MatrixWorkspace>;

namespace Mantid
{
namespace Kernel
{

template<> DLLExport
Mantid::API::MatrixWorkspace_sptr PropertyManager::getValue<Mantid::API::MatrixWorkspace_sptr>(const std::string &name) const
{
  PropertyWithValue<Mantid::API::MatrixWorkspace_sptr>* prop =
                    dynamic_cast<PropertyWithValue<Mantid::API::MatrixWorkspace_sptr>*>(getPointerToProperty(name));
  if (prop)
  {
    return *prop;
  }
  else
  {
    std::string message = "Attempt to assign property "+ name +" to incorrect type";
    throw std::runtime_error(message);
  }
}

template<> DLLExport
Mantid::API::MatrixWorkspace_const_sptr PropertyManager::getValue<Mantid::API::MatrixWorkspace_const_sptr>(const std::string &name) const
{
  PropertyWithValue<Mantid::API::MatrixWorkspace_sptr>* prop =
                    dynamic_cast<PropertyWithValue<Mantid::API::MatrixWorkspace_sptr>*>(getPointerToProperty(name));
  if (prop)
  {
    return prop->operator()();
  }
  else
  {
    std::string message = "Attempt to assign property "+ name +" to incorrect type";
    throw std::runtime_error(message);
  }
}


} // namespace Kernel
} // namespace Mantid

///\endcond TEMPLATE
