// LoadNeXus
// @author Freddie Akeroyd, STFC ISIS Faility
// @author Ronald Fowler, e_Science  - updated to be wrapper to either LoadMuonNeuxs or LoadIsisNexus
//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidNexus/LoadNeXus.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidNexus/NeXusUtils.h"
#include "MantidKernel/ArrayProperty.h"

#include <cmath>
#include <boost/shared_ptr.hpp>

namespace Mantid
{
namespace NeXus
{

  // Register the algorithm into the algorithm factory
  DECLARE_ALGORITHM(LoadNeXus)

  using namespace Kernel;
  using namespace API;
  using namespace DataObjects;

  Logger& LoadNeXus::g_log = Logger::get("LoadNeXus");

  /// Empty default constructor
  LoadNeXus::LoadNeXus() :
      Algorithm(), m_filename()
  {
  }

  /** Initialisation method.
   *
   */
  void LoadNeXus::init()
  {
    // Declare required input parameters for all sub algorithms
    std::vector<std::string> exts;
    exts.push_back("NXS");
    exts.push_back("nxs");
    declareProperty("Filename","",new FileValidator(exts));
    declareProperty(new WorkspaceProperty<Workspace2D>("OutputWorkspace","",Direction::Output));

	// Declare optional input parameters
    BoundedValidator<int> *mustBePositive = new BoundedValidator<int>();
    mustBePositive->setLower(0);
    declareProperty("spectrum_min",0, mustBePositive);
    declareProperty("spectrum_max",0, mustBePositive->clone());
  
	declareProperty(new ArrayProperty<int>("spectrum_list"));
  }

  /** Executes the algorithm. Reading in the file and creating and populating
   *  the output workspace
   *
   *  @throw runtime_error Thrown if algorithm cannot execute
   */
  void LoadNeXus::exec()
  {
    // Retrieve the filename and output workspace name from the properties
    m_filename = getPropertyValue("Filename");
    m_workspace = getPropertyValue("OutputWorkspace");

	// Test the given filename to see if it contains the field "analysis" with value "muonTD"
	// within the first NXentry.
	// If so, assume it is a Muon Nexus file (version 1) and pass to the LoadMuonNexus algorithm
	// Otherwise try LoadIsisNexus.
	std::string dataName="analysis", muonTd="muonTD";
	std::string value;
    NeXusUtils *nexusFile= new NeXusUtils();
    std::vector<std::string> entryName,definition;
    int count=nexusFile->getNexusEntryTypes(m_filename,entryName,definition);
	if( definition[0]==muonTd )
	{
		runLoadMuonNexus();
	}
	else if( entryName[0]=="mantid_workspace_1" )
	{
		runLoadNexusProcessed();
	}
	else
	{
        g_log.error("Unable to open file " + m_filename);
        throw Exception::FileError("Unable to open File:" , m_filename);
	}
    return;
  }

  void LoadNeXus::runLoadMuonNexus()
  {
      Algorithm_sptr loadMuonNexus = createSubAlgorithm("LoadMuonNexus");
      // Pass through the same input filename
      loadMuonNexus->setPropertyValue("Filename",m_filename);
      // Set the workspace property
	  std::string outputWorkspace="OutputWorkspace";
      loadMuonNexus->setPropertyValue(outputWorkspace,m_workspace);
	  //
      Property *specList = getProperty("spectrum_list");
      if( !(specList->isDefault()) )
	     loadMuonNexus->setPropertyValue("spectrum_list",getPropertyValue("spectrum_list"));
	  //
      Property *specMax = getProperty("spectrum_max");
      if( !(specMax->isDefault()) )
	  {
	     loadMuonNexus->setPropertyValue("spectrum_max",getPropertyValue("spectrum_max"));
	     loadMuonNexus->setPropertyValue("spectrum_min",getPropertyValue("spectrum_min"));
	  }

      // Now execute the sub-algorithm. Catch and log any error, but don't stop.
      try
      {
        loadMuonNexus->execute();
      }
      catch (std::runtime_error& err)
      {
        g_log.error("Unable to successfully run LoadMuonNexus sub-algorithm");
      }
      if ( ! loadMuonNexus->isExecuted() ) g_log.error("Unable to successfully run LoadMuonNexus sub-algorithm");
      // Get pointer to the workspace created
      m_localWorkspace=loadMuonNexus->getProperty(outputWorkspace); 
      setProperty<Workspace2D_sptr>(outputWorkspace,m_localWorkspace);
	  //
	  // copy pointers to any new output workspaces created by alg LoadMuonNexus to alg LoadNexus
	  // Loop through names of form "OutputWorkspace<n>" where <n> is integer from 2 upwards
	  // until name not found
	  //
	  int period=0;
	  bool noError=true;
	  while(noError)
	  {
          std::stringstream suffix;
		  period++;
          suffix << (period+1);
		  std::string opWS = outputWorkspace + suffix.str();
          std::string WSName = m_workspace + "_" + suffix.str();
		  try
		  {
              m_localWorkspace=loadMuonNexus->getProperty(opWS); 
              declareProperty(new WorkspaceProperty<DataObjects::Workspace2D>(opWS,WSName,Direction::Output));
              setProperty<Workspace2D_sptr>(opWS,m_localWorkspace);
		  }
		  catch (Exception::NotFoundError)
		  {
			  noError=false;
		  }
	  }
  }

  void LoadNeXus::runLoadNexusProcessed()
  {
  }


} // namespace NeXus
} // namespace Mantid
