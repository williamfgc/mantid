#ifndef BINARYFILE_H_
#define BINARYFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "MantidKernel/System.h"
#include "Poco/File.h"
#include "Poco/Path.h"

namespace Mantid
{
namespace Kernel
{

using std::cout;
using std::endl;
using std::ifstream;
using std::runtime_error;
using std::stringstream;
using std::string;
using std::vector;

/// Default number of items to read in from any of the files.
static const size_t DEFAULT_BLOCK_SIZE = 1000000; // 100,000

/**
 * The BinaryFile template is a helper function for loading simple binary files.
 *  - The file format must be a simple sequence of objects of type T.
 *
 */
template<typename T>
class DLLExport BinaryFile
{
public:
  //------------------------------------------------------------------------------------
  /// Empty constructor
  BinaryFile()
  { }

  //------------------------------------------------------------------------------------
  /// Constructor - open a file
  BinaryFile(std::string filename)
  {
    this->open(filename);
  }

  //------------------------------------------------------------------------------------
  /** Open a file and keep a handle to the file
   * @param filename full path to open
   * @throws runtime_error if the file size is not an even multiple of the type size
   * @throws invalid_argument if the file does not exist
   * */
  void open(std::string filename)
  {
    if (!Poco::File(filename).exists())
    {
      stringstream msg;
      msg << "BinaryFile::open: File " << filename << " was not found.";
      throw std::invalid_argument("File does not exist.");
    }
    //Open the file
    this->handle = new ifstream(filename.c_str(), std::ios::binary);
    //Count the # of elements.
    this->num_elements = this->getFileSize();
    //Make sure we are starting at 0
    this->offset = 0;

  }

  //------------------------------------------------------------------------------------
  /** Close the file
   * */
  void close()
  {
    delete handle;
  }


  //-----------------------------------------------------------------------------
  /** Get the size of a file as a multiple of a particular data type
   * @throw runtime_error if the file size is not compatible
   * @throw runtime_error if the handle is not open.
   * */
  size_t getFileSize()
  {
    this->obj_size = sizeof(T);

    if (!handle) {
      throw runtime_error("BinaryFile::getFileSize: Cannot find the size of a file from a null handle");
    }

    // get the size of the file in bytes and reset the handle back to the beginning
    handle->seekg(0, std::ios::end);
    size_t filesize = handle->tellg();
    handle->seekg(0, std::ios::beg);

    // check the file is a compatible size
    if (filesize % obj_size != 0) {
      stringstream msg;
      msg << "BinaryFile::getFileSize: File size is not compatible with data size ";
      msg << filesize << "%" << obj_size << "=";
      msg << filesize % obj_size;
      throw runtime_error(msg.str());
    }

    return filesize / sizeof(T);
  }

  //-----------------------------------------------------------------------------
  /// Returns the # of elements in the file
  size_t getNumElements()
  {
    return this->num_elements;
  }


  //-----------------------------------------------------------------------------
  /** Get a buffer size for loading blocks of data.
   * @param num_items
   */
  size_t getBufferSize(const size_t num_items)
  {
    if (num_items < DEFAULT_BLOCK_SIZE)
      return num_items;
    else
      return DEFAULT_BLOCK_SIZE;
  }

  //-----------------------------------------------------------------------------
  /**
   * Loads the entire contents of the file into a pointer to a std::vector.
   */
  std::vector<T> * loadAll()
  {
    //Initialize the pointer
    std::vector<T> * data = new std::vector<T>;

    //A buffer to load from
    size_t buffer_size = getBufferSize(num_elements);
    T * buffer = new T[buffer_size];

    //Make sure we are at the beginning
    offset = 0;
    handle->seekg(0, std::ios::beg);

    size_t loaded_size;
    while (offset < num_elements)
    {
      //Load that block of data
      loaded_size = loadBlock(buffer, buffer_size);
      // Insert into the data
      data->insert(data->end(), buffer, (buffer + loaded_size));
    }

//    while (offset < num_elements) {
//      // read a section and put it into the object
//      handle->read(reinterpret_cast<char *>(buffer), buffer_size * obj_size);
//      //Insert into the data
//      data->insert(data->end(), buffer, (buffer + buffer_size));
//      //We are further along
//      offset += buffer_size;
//
//      // make sure not to read past EOF
//      if (offset + buffer_size > num_elements)
//        buffer_size = num_elements - offset;
//    }

    //Here's your vector!
    return data;
  }


  //-----------------------------------------------------------------------------
  /**
   * Loads a single block from file and returns a pointer to a vector containing it.
   *  This can be called repeatedly to load an entire file.
   *
   * @param block_size: how many elements to load in the block. If there are not enough elements,
   *  the vector returned is smaller than block_size
   * @param buffer: array of block_size[] of T; must have been allocated before.
   * @retrun loaded_size, actually how many elements were loaded.
   */
  size_t loadBlock(T * buffer, size_t block_size)
  {
    size_t loaded_size;
    //Limit how much is loaded
    loaded_size = block_size;
    if (offset + loaded_size > num_elements)
      loaded_size = num_elements - offset;
    //Read it right into the buffer
    handle->read(reinterpret_cast<char *>(buffer), loaded_size * obj_size);
    offset += loaded_size;
    return loaded_size;
  }


private:
  /// File stream
  std::ifstream * handle;
  //Size of each object.
  size_t obj_size;
  /// Number of elements of size T in the file
  size_t num_elements;
  /// Offset into the file, if loading in blocks.
  size_t offset;

};


} //Namespace Kernel

} //Namespace Mantid

#endif /* BINARYFILE_H_ */
