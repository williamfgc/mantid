#include "MantidAPI/CatalogFactory.h"
#include "MantidAPI/CatalogManager.h"
#include "MantidAPI/CompositeCatalog.h"
#include "MantidKernel/ConfigService.h"
#include "MantidKernel/FacilityInfo.h"

namespace Mantid
{
  namespace API
  {
    CatalogManagerImpl::CatalogManagerImpl() : m_activeCatalogs() {}

    CatalogManagerImpl::~CatalogManagerImpl(){}

    /**
     * Creates a new catalog and adds it to the compositeCatalog and activeCatalog list.
     * @param facilityName :: The name of the facility to obtain the catalog name from.
     * @return A catalog for the facility specified.
     */
    ICatalog_sptr CatalogManagerImpl::create(const std::string &facilityName)
    {
      std::string className = Kernel::ConfigService::Instance().getFacility(facilityName).catalogInfo().catalogName();
      auto catalog = CatalogFactory::Instance().create(className);
      m_activeCatalogs.insert(std::make_pair(boost::lexical_cast<std::string>(rand() + 10),catalog));
      return catalog;
    }

    /**
     * Obtain a specific catalog using the sessionID, otherwise return all active catalogs.
     * @param sessionID :: The session to search for in the active catalogs list.
     * @return A specific catalog using the sessionID, otherwise returns all active catalogs
     */
    ICatalog_sptr CatalogManagerImpl::getCatalog(const std::string &sessionID)
    {
      if(sessionID.empty())
      {
        auto composite = boost::make_shared<CompositeCatalog>();
        for (auto item = m_activeCatalogs.begin(); item != m_activeCatalogs.end(); ++item)
        {
          composite->add(item->second);
        }
        return composite;
      }

      for(auto iter = m_activeCatalogs.begin(); iter != m_activeCatalogs.end(); ++iter)
      {
        if (sessionID == iter->first->getSessionId()) return iter->second;
      }

      // If we reached this point then the session is corrupt/invalid.
      throw std::runtime_error("The session ID you have provided is invalid");
    }

    /**
     * Destroy and remove a specific catalog from the active catalogs list and the composite catalog.
     * @param sessionID :: The session to search for in the active catalogs list.
     */
    void CatalogManagerImpl::destroyCatalog(const std::string& sessionID)
    {
      for(auto iter = m_activeCatalogs.begin(); iter != m_activeCatalogs.end(); ++iter)
      {
        if (sessionID == iter->first->getSessionId())
        {
          iter->second->logout();
          m_activeCatalogs.erase(iter);
        }
      }
    }

    /**
     * Destroy all active catalogs.
     */
    void CatalogManagerImpl::destroyCatalogs()
    {
      for(auto item = m_activeCatalogs.begin(); item != m_activeCatalogs.end(); ++item)
      {
        item->second->logout();
      }
      m_activeCatalogs.clear();
    }

  }
}
