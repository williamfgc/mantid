// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidKernel/NexusDescriptor.h"

// clang-format off
#include <nexus/NeXusFile.hpp>
#include <nexus/NeXusException.hpp>
// clang-format on

#include <Poco/File.h>
#include <Poco/Path.h>

#include <cstring>
#include <string>

#include <hdf5.h>

namespace Mantid {
namespace Kernel {
//---------------------------------------------------------------------------------------------------------------------------
// static NexusDescriptor constants
//---------------------------------------------------------------------------------------------------------------------------
/// Size of HDF magic number
const size_t NexusDescriptor::HDFMagicSize = 4;
/// HDF cookie that is stored in the first 4 bytes of the file.
const unsigned char NexusDescriptor::HDFMagic[4] = {
    '\016', '\003', '\023', '\001'}; // From HDF4::hfile.h

/// Size of HDF5 signature
size_t NexusDescriptor::HDF5SignatureSize = 8;
/// signature identifying a HDF5 file.
const unsigned char NexusDescriptor::HDF5Signature[8] = {
    137, 'H', 'D', 'F', '\r', '\n', '\032', '\n'};

namespace {
//---------------------------------------------------------------------------------------------------------------------------
// Anonymous helper methods to use isHDF methods to use an open file handle
//---------------------------------------------------------------------------------------------------------------------------

/**
 * Currently simply checks for the HDF signatures and returns true if one of
 * them is found
 * @param fileHandle A file handled opened and pointing at the start of the
 * file. On return the
 * fileHandle is left at the start of the file
 * @param version One of the NexusDescriptor::Version enumerations specifying
 * the required version
 * @return True if the file is considered hierarchical, false otherwise
 */
bool isHDFHandle(FILE *fileHandle, NexusDescriptor::Version version) {
  if (!fileHandle)
    throw std::invalid_argument(
        "HierarchicalFileDescriptor::isHierarchical - Invalid file handle");

  bool result(false);

  // HDF4 check requires 4 bytes,  HDF5 check requires 8 bytes
  // Use same buffer and waste a few bytes if only checking HDF4
  unsigned char buffer[8] = {'0', '0', '0', '0', '0', '0', '0', '0'};
  if (NexusDescriptor::HDF5SignatureSize !=
      std::fread(static_cast<void *>(&buffer), sizeof(unsigned char),
                 NexusDescriptor::HDF5SignatureSize, fileHandle)) {
    throw std::runtime_error("Error while reading file");
  }

  // Number of bytes read doesn't matter as if it is not enough then the memory
  // simply won't match
  // as the buffer has been "zeroed"
  if (version == NexusDescriptor::Version5 ||
      version == NexusDescriptor::AnyVersion) {
    result = (std::memcmp(&buffer, &NexusDescriptor::HDF5Signature,
                          NexusDescriptor::HDF5SignatureSize) == 0);
  }
  if (!result && (version == NexusDescriptor::Version4 ||
                  version == NexusDescriptor::AnyVersion)) {
    result = (std::memcmp(&buffer, &NexusDescriptor::HDFMagic,
                          NexusDescriptor::HDFMagicSize) == 0);
  }

  // Return file stream to start of file
  std::rewind(fileHandle);
  return result;
}
} // namespace

//---------------------------------------------------------------------------------------------------------------------------
// static NexusDescriptor methods
//---------------------------------------------------------------------------------------------------------------------------

/**
 * Checks for the HDF signatures and returns true if one of them is found
 * @param filename A string filename to check
 * @param version One of the NexusDescriptor::Version enumerations specifying
 * the required version
 * @return True if the file is considered hierarchical, false otherwise
 */
bool NexusDescriptor::isHDF(const std::string &filename,
                            const Version version) {
  FILE *fd = fopen(filename.c_str(), "rb");
  if (!fd) {
    throw std::invalid_argument(
        "HierarchicalFileDescriptor::isHierarchical - Unable to open file '" +
        filename + "'");
  }
  const bool result = isHDFHandle(fd, version); // use anonymous helper
  fclose(fd);
  return result;
}

//---------------------------------------------------------------------------------------------------------------------------
// NexusDescriptor public methods
//---------------------------------------------------------------------------------------------------------------------------
/**
 * Constructs the wrapper
 * @param filename A string pointing to an existing file
 * @throws std::invalid_argument if the file is not identified to be
 * hierarchical. This currently
 * involves simply checking for the signature if a HDF file at the start of the
 * file
 */
NexusDescriptor::NexusDescriptor(const std::string &filename)
    : m_filename(), m_extension(), m_firstEntryNameType(), m_rootAttrs(),
      m_pathsToTypes(), m_file(nullptr) {
  if (filename.empty()) {
    throw std::invalid_argument("NexusDescriptor() - Empty filename '" +
                                filename + "'");
  }
  if (!Poco::File(filename).exists()) {
    throw std::invalid_argument("NexusDescriptor() - File '" + filename +
                                "' does not exist");
  }
  try {
    initialize(filename);
  } catch (::NeXus::Exception &e) {
    throw std::invalid_argument(
        "NexusDescriptor::initialize - File '" + filename +
        "' does not look like a HDF file.\n Error was: " + e.what());
  }
}

/**
 */
NexusDescriptor::~NexusDescriptor() {}

/// Returns the name & type of the first entry in the file
const std::pair<std::string, std::string> &
NexusDescriptor::firstEntryNameType() const {
  return m_firstEntryNameType;
}

/**
 * @param name The name of an attribute
 * @return True if the attribute exists, false otherwise
 */
bool NexusDescriptor::hasRootAttr(const std::string &name) const {
  return (m_rootAttrs.count(name) == 1);
}

/**
 * @param path A string giving a path using UNIX-style path separators (/), e.g.
 * /raw_data_1, /entry/bank1
 * @return True if the path exists in the file, false otherwise
 */
bool NexusDescriptor::pathExists(const std::string &path) const {
  return (m_pathsToTypes.find(path) != m_pathsToTypes.end());
}

/**
 * @param path A string giving a path using UNIX-style path separators (/), e.g.
 * /raw_data_1, /entry/bank1
 * @param type A string specifying the required type
 * @return True if the path exists in the file, false otherwise
 */
bool NexusDescriptor::pathOfTypeExists(const std::string &path,
                                       const std::string &type) const {
  auto it = m_pathsToTypes.find(path);
  if (it != m_pathsToTypes.end()) {
    return (it->second == type);
  } else
    return false;
}

/**
 * @param type A string specifying the required type
 * @return path A string giving a path using UNIX-style path separators (/),
 * e.g. /raw_data_1, /entry/bank1
 */
std::string NexusDescriptor::pathOfType(const std::string &type) const {
  auto iend = m_pathsToTypes.end();
  for (auto it = m_pathsToTypes.begin(); it != iend; ++it) {
    if (type == it->second)
      return it->first;
  }
  return "";
}

/**
 * @param classType A string name giving a class type
 * @return True if the type exists in the file, false otherwise
 */
bool NexusDescriptor::classTypeExists(const std::string &classType) const {
  auto iend = m_pathsToTypes.end();
  for (auto it = m_pathsToTypes.begin(); it != iend; ++it) {
    if (classType == it->second)
      return true;
  }
  return false;
}

//---------------------------------------------------------------------------------------------------------------------------
// NexusDescriptor private methods
//---------------------------------------------------------------------------------------------------------------------------

/**
 * Creates the internal cached structure of the file as a tree of nodes
 */
void NexusDescriptor::initialize(const std::string &filename) {
  m_filename = filename;
  m_extension = "." + Poco::Path(filename).getExtension();

  m_file = std::make_unique<::NeXus::File>(this->filename());

  m_file->openPath("/");
  m_rootAttrs.clear();
  m_pathsToTypes.clear();
  walkFile(*m_file, "", "", m_pathsToTypes, 0);
}

/**
 * Cache the structure in the given maps
 * @param file An open NeXus File object
 * @param rootPath The current path that is open in the file
 * @param className The class of the current open path
 * @param pmap [Out] An output map filled with mappings of path->type
 * @param level An integer defining the current level in the file
 */
void NexusDescriptor::walkFile(::NeXus::File &file, const std::string &rootPath,
                               const std::string &className,
                               std::map<std::string, std::string> &pmap,
                               int level) {
  if (!rootPath.empty()) {
    pmap.emplace(rootPath, className);
  }
  if (level == 0) {
    auto attrInfos = file.getAttrInfos();
    for (auto &attrInfo : attrInfos) {
      m_rootAttrs.insert(attrInfo.name);
    }
  }

  auto dirents = file.getEntries();
  auto itend = dirents.end();
  for (auto it = dirents.begin(); it != itend; ++it) {
    const std::string &entryName = it->first;
    const std::string &entryClass = it->second;
    const std::string entryPath =
        std::string(rootPath).append("/").append(entryName);
    if (entryClass == "SDS" || entryClass == "ILL_data_scan_vars") {
      pmap.emplace(entryPath, entryClass);
    } else if (entryClass == "CDF0.0") {
      // Do nothing with this
    } else {
      if (level == 0)
        m_firstEntryNameType = (*it); // copy first entry name & type
      file.openGroup(entryName, entryClass);
      walkFile(file, entryPath, entryClass, pmap, level + 1);
    }
  }
  file.closeGroup();
}

namespace {

herr_t readStringAttribute(hid_t attr, char **data) {
  herr_t iRet = 0;
  hid_t atype = -1;
  hid_t space;
  int ndims;
  hsize_t thedims[H5S_MAX_RANK], sdim;

  atype = H5Aget_type(attr);
  sdim = H5Tget_size(atype);
  space = H5Aget_space(attr);
  ndims = H5Sget_simple_extent_dims(space, thedims, NULL);

  if (ndims == 0) {
    if (H5Tis_variable_str(atype)) {
      hid_t btype = H5Tget_native_type(atype, H5T_DIR_ASCEND);
      iRet = H5Aread(attr, btype, data);
      H5Tclose(btype);
    } else {
      *data = (char *)malloc(sdim + 1);
      iRet = H5Aread(attr, atype, *data);
      (*data)[sdim] = '\0';
    }
  } else if (ndims == 1) {
    unsigned int i;
    char **strings;

    strings = (char **)malloc(thedims[0] * sizeof(char *));

    if (!H5Tis_variable_str(atype)) {
      strings[0] = (char *)malloc(thedims[0] * sdim * sizeof(char));
      for (i = 1; i < thedims[0]; i++) {
        strings[i] = strings[0] + i * sdim;
      }
    }

    iRet = H5Aread(attr, atype, strings[0]);
    *data = (char *)calloc((sdim + 2) * thedims[0], sizeof(char));
    for (i = 0; i < thedims[0]; i++) {
      if (i == 0) {
        strncpy(*data, strings[i], sdim);
      } else {
        strcat(*data, ", ");
        strncat(*data, strings[i], sdim);
      }
    }
    if (H5Tis_variable_str(atype)) {
      H5Dvlen_reclaim(atype, space, H5P_DEFAULT, strings);
    } else {
      free(strings[0]);
    }

    free(strings);
  } else {
    *data = strdup(" higher dimensional string array");
  }

  H5Tclose(atype);
  H5Sclose(space);
  if (iRet < 0)
    return NX_ERROR;
  return NX_OK;
}

herr_t readStringAttributeN(hid_t attr, char *data, int maxlen) {
  herr_t iRet;
  char *vdat = NULL;
  iRet = readStringAttribute(attr, &vdat);
  if (iRet >= 0) {
    strncpy(data, vdat, maxlen);
    free(vdat);
  }
  data[maxlen - 1] = '\0';
  return iRet;
}

void getGroup(hid_t groupID,
              std::map<std::string, std::set<std::string>> &allEntries) {

  auto lf_getNxClassAttribute = [&](hid_t groupID,
                                    const char *objectName) -> std::string {
    std::string attribute = "";
    hid_t attributeID = H5Aopen_by_name(groupID, objectName, "NX_class",
                                        H5P_DEFAULT, H5P_DEFAULT);
    if (attributeID < 0) {
      H5Aclose(attributeID);
      return attribute;
    }

    hid_t type = H5T_C_S1;
    hid_t atype = H5Tcopy(type);
    char data[128];
    H5Tset_size(atype, sizeof(data));
    readStringAttributeN(attributeID, data, sizeof(data));
    attribute = std::string(data);
    H5Tclose(atype);
    H5Aclose(attributeID);

    return attribute;
  };

  // using HDF5 C API
  constexpr size_t maxLength = 1024;
  char groupName[maxLength];
  char memberName[maxLength];
  ssize_t groupNameLength = H5Iget_name(groupID, groupName, maxLength);
  hsize_t nObjects = 0;
  // FIXME handler errors
  H5Gget_num_objs(groupID, &nObjects);

  const std::string groupNameStr(groupName, groupNameLength);
  const std::string nxClass =
      (groupNameStr == "/")
          ? ""
          : lf_getNxClassAttribute(groupID, groupNameStr.c_str());

  // std::cout << " GROUP: " << groupName << " nObjects: " << nObjects << "\n";

  if (!nxClass.empty()) {
    allEntries[nxClass].insert(groupNameStr);
    // std::cout << " ENTRY " << nxClass << " " << groupNameStr << "\n";
  }

  for (unsigned int i = 0; i < nObjects; ++i) {

    const int type = H5Gget_objtype_by_idx(groupID, static_cast<size_t>(i));
    ssize_t memberNameLength = H5Gget_objname_by_idx(
        groupID, static_cast<hsize_t>(i), memberName, maxLength);

    if (type == H5O_TYPE_GROUP) {
      hid_t subGroupID = H5Gopen2(groupID, memberName, H5P_DEFAULT);
      getGroup(subGroupID, allEntries);
      // FIXME handle errors
      H5Gclose(subGroupID);

    } else if (type == H5O_TYPE_DATASET) {
      const std::string memberNameStr(memberName, memberNameLength);
      const std::string absoluteEntryName = groupNameStr + "/" + memberNameStr;
      allEntries["SDS"].insert(absoluteEntryName);
      // std::cout << " ENTRY SDS " << absoluteEntryName << "\n";
    }
  }
}

} // namespace

void NexusDescriptor::setAllEntries() {

  hid_t fileID = H5Fopen(m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (fileID < 0) {
    throw std::invalid_argument("ERROR: couldn't open hdf5 file " + m_filename +
                                "\n");
  }

  // FIXME here add a check to file to verify it's a HDF5 file
  hid_t groupID = H5Gopen2(fileID, "/", H5P_DEFAULT);
  getGroup(groupID, m_allEntries);
  // FIXME handle errors
  H5Gclose(groupID);
  H5Fclose(fileID);
}

} // namespace Kernel
} // namespace Mantid
