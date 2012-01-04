//---------------------------------------
// Includes
//---------------------------------------
#include "MantidPythonAPI/WorkspaceProxies.h"
#include "MantidPythonAPI/MantidVecHelper.h"

#include "MantidAPI/WorkspaceOpOverloads.h"
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/WorkspaceGroup.h"
#include "MantidAPI/IMDWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/IMDHistoWorkspace.h"

using namespace Mantid::API;

namespace Mantid
{
  namespace PythonAPI
  {
    using API::MatrixWorkspace_sptr;
    using API::WorkspaceGroup_sptr;

    //****************************************
    //
    // MatrixWorkspaceWrapper
    //
    //****************************************
    /**
    * Returns the X values from the spectra at the given index wrapped in a read-only numpy array
    * @param self :: A reference to the workspace
    * @param index :: The index of the workspace
    * @returns A numpy array for the given index
    */
    PyObject * MatrixWorkspaceWrapper::readX(API::MatrixWorkspace& self, int index)
    {
      return MantidVecHelper::createPythonWrapper(self.readX(index), true);
    }

    /**
    * Returns the Y values from the spectra at the given index wrapped in a read-only numpy array
    * @param self :: A reference to the workspace
    * @param index :: The index of the workspace
    * @returns A numpy array for the given index
    */
    PyObject * MatrixWorkspaceWrapper::readY(API::MatrixWorkspace& self, int index)
    {
      return MantidVecHelper::createPythonWrapper(self.readY(index), true);
    }

    /**
    * Returns the E values from the spectra at the given index wrapped in a read-only numpy array
    * @param self :: A reference to the workspace
    * @param index :: The index of the workspace
    * @returns A numpy array for the given index
    */
    PyObject * MatrixWorkspaceWrapper::readE(API::MatrixWorkspace& self, int index)
    {
      return MantidVecHelper::createPythonWrapper(self.readE(index), true);
    }

    /**
    * Returns the Dx values from the spectra at the given index wrapped in a read-only numpy array
    * @param self :: A reference to the workspace
    * @param index :: The index of the workspace
    * @returns A numpy array for the given index
    */
    PyObject * MatrixWorkspaceWrapper::readDx(API::MatrixWorkspace& self, int index)
    {
      return MantidVecHelper::createPythonWrapper(self.readDx(index), true);
    }

    //*********************************************************************************
    //
    // Binary Operation helpers
    //
    //**********************************************************************************

    /** Binary operation for two workspaces. Generic for IMDWorkspaces or MatrixWorkspaces...
     * Called by python overloads for _binary_op (see api_exports.cpp)
     *
     * @param lhs :: the left hand side workspace of the operation
     * @param rhs :: the right hand side workspace of the operation
     * @param op :: The operation
     * @param name :: The output name
     * @param inplace :: is this is an inplace operation (i.e. does the output overwrite the lhs
     * @param reverse :: Unused parameter. Here for consistent interface
     * @returns The resulting workspace
     */
    template<typename LHSType, typename RHSType, typename ResultType>
    ResultType performBinaryOp(const LHSType lhs, const RHSType rhs,
                        const std::string& op, const std::string & name,
                        bool inplace,bool reverse)
    {
      std::string algoName = op;

      // ----- Determine which version of the algo should be called -----
      MatrixWorkspace_const_sptr lhs_mat = boost::dynamic_pointer_cast<const MatrixWorkspace>(lhs);
      MatrixWorkspace_const_sptr rhs_mat = boost::dynamic_pointer_cast<const MatrixWorkspace>(rhs);
      WorkspaceGroup_const_sptr lhs_grp = boost::dynamic_pointer_cast<const WorkspaceGroup>(lhs);
      WorkspaceGroup_const_sptr rhs_grp = boost::dynamic_pointer_cast<const WorkspaceGroup>(rhs);

      if ( (lhs_mat || lhs_grp) && (rhs_mat || rhs_grp) )
        // Both sides are matrixworkspace - use the original algos (e..g "Plus.")
        algoName = op;
      else
        // One of the workspaces must be MDHistoWorkspace or MDEventWorkspace
        // Use the MD version, e.g. "PlusMD"
        algoName = op + "MD";

      ResultType result;
      std::string error("");
      try
      {
        if( reverse )
        {
          result = API::OperatorOverloads::executeBinaryOperation<RHSType, LHSType, ResultType>(algoName, rhs, lhs, inplace, false, name, true);
        }
        else
        {
          result = API::OperatorOverloads::executeBinaryOperation<LHSType, RHSType, ResultType>(algoName, lhs, rhs, inplace, false, name, true);
        }
      }
      catch(std::runtime_error & exc)
      {
        error = exc.what();
        if( error.find("algorithm") == 0 )
        {
          error = "Unknown binary operation requested: " + op;
          throw std::runtime_error(error);
        }
        else
        {
          throw;
        }
      }
      return result;
    }


    /** 
    * Perform the given binary operation on a workspace and a double.
    * Generic to MDWorkspaces.
    * Called by python overloads for _binary_op (see api_exports.cpp)
    *
    * @param inputWS :: The input workspace
    * @param value :: The input value
    * @param op :: The operation
    * @param name :: The output name
    * @param inplace :: If true, then the lhs argument is replaced by the result of the operation.
    * @param reverse :: If true then the double is the lhs argument
    * @return A shared pointer to the result workspace
    */
    template<typename LHSType, typename ResultType>
    ResultType performBinaryOpWithDouble(const LHSType inputWS, const double value,
        const std::string& op, const std::string & name,
        bool inplace, bool reverse)
    {
      std::string algoName = op;

      // Create the single valued workspace first so that it is run as a top-level algorithm
      // such that it's history can be recreated
      API::Algorithm_sptr alg = API::AlgorithmManager::Instance().createUnmanaged("CreateSingleValuedWorkspace");
      alg->setChild(false);
      alg->initialize();
      alg->setProperty<double>("DataValue",value);
      const std::string & tmp_name("__tmp_binary_operation_double");
      alg->setPropertyValue("OutputWorkspace", tmp_name);
      alg->execute();
      MatrixWorkspace_sptr singleValue;
      API::AnalysisDataServiceImpl & data_store = API::AnalysisDataService::Instance();
      if( alg->isExecuted() )
      {
        singleValue = boost::dynamic_pointer_cast<API::MatrixWorkspace>(data_store.retrieve(tmp_name));
      }
      else
      {
        throw std::runtime_error("performBinaryOp: Error in execution of CreateSingleValuedWorkspace");
      }      
      // Call the function above with the signle-value workspace
      ResultType result = performBinaryOp<LHSType, MatrixWorkspace_sptr, ResultType>(inputWS, singleValue, algoName, name, inplace, reverse);
      // Delete the temporary
      data_store.remove(tmp_name);
      return result;
    }



    // Concrete instantations
    template IMDWorkspace_sptr performBinaryOp(const IMDWorkspace_sptr, const IMDWorkspace_sptr, const std::string& , const std::string & name,
        bool, bool);
    template MatrixWorkspace_sptr performBinaryOp(const MatrixWorkspace_sptr, const MatrixWorkspace_sptr, const std::string& , const std::string & name,
        bool, bool);
    template WorkspaceGroup_sptr performBinaryOp(const IMDWorkspace_sptr, const WorkspaceGroup_sptr, const std::string& , const std::string & name,
        bool, bool);
    template WorkspaceGroup_sptr performBinaryOp(const WorkspaceGroup_sptr, const IMDWorkspace_sptr, const std::string& , const std::string & name,
        bool, bool);
    template WorkspaceGroup_sptr performBinaryOp(const WorkspaceGroup_sptr, const WorkspaceGroup_sptr, const std::string& , const std::string & name,
        bool, bool);

    template IMDHistoWorkspace_sptr performBinaryOp(const IMDHistoWorkspace_sptr, const IMDHistoWorkspace_sptr, const std::string& , const std::string & name,
        bool, bool);
    template IMDHistoWorkspace_sptr performBinaryOp(const IMDHistoWorkspace_sptr, const MatrixWorkspace_sptr, const std::string& , const std::string & name,
        bool, bool);


    // Double variants
    template IMDWorkspace_sptr performBinaryOpWithDouble(const IMDWorkspace_sptr, const double, const std::string& op,
        const std::string &, bool, bool);
    template MatrixWorkspace_sptr performBinaryOpWithDouble(const MatrixWorkspace_sptr, const double, const std::string& op,
        const std::string &, bool, bool);
    template IMDHistoWorkspace_sptr performBinaryOpWithDouble(const IMDHistoWorkspace_sptr, const double, const std::string& op,
        const std::string &, bool, bool);
    template WorkspaceGroup_sptr performBinaryOpWithDouble(const WorkspaceGroup_sptr, const double, const std::string& op,
        const std::string &, bool, bool);

  }

}
