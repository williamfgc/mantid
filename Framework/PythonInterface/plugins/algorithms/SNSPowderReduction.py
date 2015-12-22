#pylint: disable=invalid-name,no-init,too-many-lines
import os
import numpy

import mantid
import mantid.api
import mantid.simpleapi as api
from mantid.api import *
from mantid.kernel import *

if AlgorithmFactory.exists('GatherWorkspaces'):
    HAVE_MPI = True
    from mpi4py import MPI
    mpiRank = MPI.COMM_WORLD.Get_rank()
else:
    HAVE_MPI = False
    mpiRank = 0 # simplify if clauses


EVENT_WORKSPACE_ID = "EventWorkspace"

def noRunSpecified(runs):
    """

    :param runs:
    :return:
    """
    # TODO/NOW - Doc!
    assert(runs, numpy.ndarray)
    if runs.size <= 0:
        return True
    if runs.size == 1:
        return runs[0] <= 0
    return False

def allEventWorkspaces(*args):
    result = True

    for arg in args:
        assert isinstance(arg, str)
        workspace = AnalysisDataService.retrieve(arg)
        result = result and (workspace.id() == EVENT_WORKSPACE_ID)

    return result

#pylint: disable=too-many-instance-attributes
class SNSPowderReduction(DataProcessorAlgorithm):
    COMPRESS_TOL_TOF = .01
    _resampleX = None
    _binning = None
    _bin_in_dspace = None
    _instrument = None
    _filterBadPulses = None
    _removePromptPulseWidth = None
    _LRef = None
    _DIFCref = None
    _wavelengthMin = None
    _wavelengthMax = None
    _vanPeakFWHM = None
    _vanSmoothing = None
    _vanRadius = None
    _scaleFactor = None
    _outDir = None
    _outPrefix = None
    _outTypes = None
    _chunks = None
    _splitws = None
    _splitinfotablews = None
    _normalisebycurrent = None
    _lowResTOFoffset = None
    _focusPos = None
    _charTable = None
    iparmFile = None
    _info = None

    def category(self):
        return "Diffraction\\Reduction"

    def name(self):
        return "SNSPowderReduction"

    def summary(self):
        " "
        return "The algorithm used for reduction of powder diffraction data obtained on SNS instruments (e.g. PG3) "

    def PyInit(self):
        sns = ConfigService.getFacility("SNS")
        instruments = []
        for item in sns.instruments("Neutron Diffraction"):
            instruments.append(item.shortName())
        self.declareProperty("Instrument", "PG3", StringListValidator(instruments), "Powder diffractometer's name")
        arrvalidator = IntArrayBoundedValidator()
        arrvalidator.setLower(0)
        self.declareProperty(IntArrayProperty("RunNumber", values=[0], validator=arrvalidator,\
                             direction=Direction.Input), "Number of sample run or 0 for only Vanadium and/or Background")
        extensions = [ "_histo.nxs", "_event.nxs", "_runinfo.xml", ".nxs.h5"]
        self.declareProperty("Extension", "_event.nxs",
                             StringListValidator(extensions))
        self.declareProperty("PreserveEvents", True,
                             "Argument to supply to algorithms that can change from events to histograms.")
        self.declareProperty("Sum", False,
                             "Sum the runs. Does nothing for characterization runs")
        self.declareProperty("PushDataPositive", "None",
                             StringListValidator(["None", "ResetToZero", "AddMinimum"]),
                             "Add a constant to the data that makes it positive over the whole range.")
        arrvalidatorBack = IntArrayBoundedValidator()
        arrvalidator.setLower(-1)
        self.declareProperty(IntArrayProperty("BackgroundNumber", values=[0], validator=arrvalidatorBack),
                             doc="If specified overrides value in CharacterizationRunsFile If -1 turns off correction.")
        arrvalidatorVan = IntArrayBoundedValidator()
        arrvalidator.setLower(-1)
        self.declareProperty(IntArrayProperty("VanadiumNumber", values=[0], validator=arrvalidatorVan),
                             doc="If specified overrides value in CharacterizationRunsFile. If -1 turns off correction.")
        arrvalidatorVanBack = IntArrayBoundedValidator()
        arrvalidator.setLower(-1)
        self.declareProperty(IntArrayProperty("VanadiumBackgroundNumber", values=[0], validator=arrvalidatorVanBack),
                             doc="If specified overrides value in CharacterizationRunsFile. If -1 turns off correction.")
        self.declareProperty(FileProperty(name="CalibrationFile",defaultValue="",action=FileAction.Load,\
                                      extensions = [".h5", ".hd5", ".hdf", ".cal"]))
        self.declareProperty(FileProperty(name="CharacterizationRunsFile",defaultValue="",action=FileAction.OptionalLoad,\
                                      extensions = ["txt"]),"File with characterization runs denoted")
        self.declareProperty(FileProperty(name="ExpIniFilename", defaultValue="", action=FileAction.OptionalLoad,
                                          extensions=[".ini"]))
        self.declareProperty("UnwrapRef", 0.,
                             "Reference total flight path for frame unwrapping. Zero skips the correction")
        self.declareProperty("LowResRef", 0.,
                             "Reference DIFC for resolution removal. Zero skips the correction")
        self.declareProperty("CropWavelengthMin", 0.,
                             "Crop the data at this minimum wavelength. Overrides LowResRef.")
        self.declareProperty("CropWavelengthMax", 0.,
                             "Crop the data at this maximum wavelength. Forces use of CropWavelengthMin.")
        self.declareProperty("RemovePromptPulseWidth", 0.0,
                             "Width of events (in microseconds) near the prompt pulse to remove. 0 disables")
        self.declareProperty("MaxChunkSize", 0.0, "Specify maximum Gbytes of file to read in one chunk.  Default is whole file.")
        self.declareProperty("FilterCharacterizations", False,
                             "Filter the characterization runs using above parameters. This only works for event files.")
        self.declareProperty(FloatArrayProperty("Binning", values=[0.,0.,0.],\
                             direction=Direction.Input), "Positive is linear bins, negative is logorithmic")
        self.declareProperty("ResampleX", 0,
                             "Number of bins in x-axis. Non-zero value overrides \"Params\" property. "+\
                             "Negative value means logorithmic binning.")
        self.declareProperty("BinInDspace", True,
                             "If all three bin parameters a specified, whether they are in dspace (true) or time-of-flight (false)")
        self.declareProperty("StripVanadiumPeaks", True,
                             "Subtract fitted vanadium peaks from the known positions.")
        self.declareProperty("VanadiumFWHM", 7, "Default=7")
        self.declareProperty("VanadiumPeakTol", 0.05,
                             "How far from the ideal position a vanadium peak can be during StripVanadiumPeaks. "\
                             "Default=0.05, negative turns off")
        self.declareProperty("VanadiumSmoothParams", "20,2", "Default=20,2")
        self.declareProperty("VanadiumRadius", .3175, "Radius for MultipleScatteringCylinderAbsorption")
        self.declareProperty("BackgroundSmoothParams", "", "Default=off, suggested 20,2")
        self.declareProperty("FilterBadPulses", 95.,
                             doc="Filter out events measured while proton charge is more than 5% below average")
        self.declareProperty("ScaleData", defaultValue=1., validator=FloatBoundedValidator(lower=0., exclusive=True),
                             doc="Constant to multiply the data before writing out. This does not apply to PDFgetN files.")
        self.declareProperty("SaveAs", "gsas",
                             "List of all output file types. Allowed values are 'fullprof', 'gsas', 'nexus', 'pdfgetn', and 'topas'")
        self.declareProperty("OutputFilePrefix", "", "Overrides the default filename for the output file (Optional).")
        self.declareProperty(FileProperty(name="OutputDirectory",defaultValue="",action=FileAction.Directory))
        self.declareProperty("FinalDataUnits", "dSpacing", StringListValidator(["dSpacing","MomentumTransfer"]))

        tableprop = ITableWorkspaceProperty("SplittersWorkspace", "", Direction.Input, PropertyMode.Optional)
        self.declareProperty(tableprop, "Splitters workspace for split event workspace.")
        infotableprop = ITableWorkspaceProperty("SplitInformationWorkspace", "", Direction.Input, PropertyMode.Optional)
        self.declareProperty(infotableprop, "Name of table workspace containing information for splitters.")

        self.declareProperty("LowResolutionSpectraOffset", -1,
                             "If larger and equal to 0, then process low resolution TOF and offset is the spectra number. "+\
                             "Otherwise, ignored.")

        self.declareProperty("NormalizeByCurrent", True, "Normalize by current")

        self.declareProperty("CompressTOFTolerance", 0.01, "Tolerance to compress events in TOF.")

        self.declareProperty(StringArrayProperty("FrequencyLogNames", ["SpeedRequest1", "Speed1", "frequency"],\
            direction=Direction.Input),\
            "Possible log names for frequency.")

        self.declareProperty(StringArrayProperty("WaveLengthLogNames", ["LambdaRequest", "lambda"],\
            direction=Direction.Input),\
            "Candidate log names for wave length.")

        return

    #pylint: disable=too-many-locals,too-many-branches,too-many-statements
    def PyExec(self):
        """ Main execution body
        """
        # get generic information
        SUFFIX = self.getProperty("Extension").value
        self._loadCharacterizations()
        self._resampleX = self.getProperty("ResampleX").value
        if self._resampleX != 0.:
            self._binning = [0.]
        else:
            self._binning = self.getProperty("Binning").value
            if len(self._binning) != 1 and len(self._binning) != 3:
                raise RuntimeError("Can only specify (width) or (start,width,stop) for binning. Found %d values." % len(self._binning))
            if len(self._binning) == 3:
                if self._binning[0] == 0. and self._binning[1] == 0. and self._binning[2] == 0.:
                    raise RuntimeError("Failed to specify the binning")
        self._bin_in_dspace = self.getProperty("BinInDspace").value
        self._instrument = self.getProperty("Instrument").value
        config['default.facility'] = "SNS"
        config['default.instrument'] = self._instrument
        self._filterBadPulses = self.getProperty("FilterBadPulses").value
        self._removePromptPulseWidth = self.getProperty("RemovePromptPulseWidth").value
        self._LRef = self.getProperty("UnwrapRef").value
        self._DIFCref = self.getProperty("LowResRef").value
        self._wavelengthMin = self.getProperty("CropWavelengthMin").value
        self._wavelengthMax = self.getProperty("CropWavelengthMax").value
        self._vanPeakFWHM = self.getProperty("VanadiumFWHM").value
        self._vanSmoothing = self.getProperty("VanadiumSmoothParams").value
        self._vanRadius = self.getProperty("VanadiumRadius").value
        calib = self.getProperty("CalibrationFile").value
        self._scaleFactor = self.getProperty("ScaleData").value
        self._outDir = self.getProperty("OutputDirectory").value
        self._outPrefix = self.getProperty("OutputFilePrefix").value.strip()
        self._outTypes = self.getProperty("SaveAs").value.lower()
        samRuns = self.getProperty("RunNumber").value
        preserveEvents = self.getProperty("PreserveEvents").value
        if HAVE_MPI and preserveEvents == True:
            self.log().warning("preserveEvents set to False for MPI tasks.")
            preserveEvents = False
        self._info = None
        self._chunks = self.getProperty("MaxChunkSize").value

        self._splitws = self.getProperty("SplittersWorkspace").value
        if self._splitws is not None:
            self.log().information("SplittersWorkspace is %s" % (str(self._splitws)))
            if len(samRuns) != 1:
                raise NotImplementedError("Reducing data with splitting cannot happen when there are more than 1 sample run.")
            timeFilterWall = self._getTimeFilterWall(self._splitws, samRuns[0], SUFFIX)
            self.log().information("The time filter wall is %s" %(str(timeFilterWall)))
        else:
            timeFilterWall = (0.0, 0.0)
            self.log().information("SplittersWorkspace is None, and thus there is NO time filter wall. ")

        self._splitinfotablews = self.getProperty("SplitInformationWorkspace").value

        self._normalisebycurrent = self.getProperty("NormalizeByCurrent").value

        # Tolerance for compress TOF event
        self.COMPRESS_TOL_TOF = float(self.getProperty("CompressTOFTolerance").value)
        if self.COMPRESS_TOL_TOF < 0.:
            self.COMPRESS_TOL_TOF = 0.01

        # Process data
        workspacelist = [] # all data workspaces that will be converted to d-spacing in the end
        samwksplist = []

        self._lowResTOFoffset = self.getProperty("LowResolutionSpectraOffset").value
        focuspos = self._focusPos
        if self._lowResTOFoffset >= 0:
            # Dealing with the parameters for editing instrument parameters
            if focuspos.has_key("PrimaryFlightPath") is True:
                l1 = focuspos["PrimaryFlightPath"]
                if l1 > 0:
                    specids = focuspos['SpectrumIDs'][:]
                    l2s = focuspos['L2'][:]
                    polars = focuspos['Polar'][:]
                    phis = focuspos['Azimuthal'][:]

                    specids.extend(specids)
                    l2s.extend(l2s)
                    polars.extend(polars)
                    phis.extend(phis)

                    focuspos['SpectrumIDs'] = specids
                    focuspos['L2'] = l2s
                    focuspos['Polar'] = polars
                    focuspos['Azimuthal'] = phis
        # ENDIF

        if self.getProperty("Sum").value:
            # Sum input sample runs and then do reduction
            if self._splitws is not None:
                raise NotImplementedError("Summing spectra and filtering events are not supported simultaneously.")

            samRun = self._focusAndSum(samRuns, SUFFIX, timeFilterWall, calib,\
                                       preserveEvents=preserveEvents)

            samRuns = [samRun]
            workspacelist.append(str(samRun))
            samwksplist.append(str(samRun))
        # ENDIF (SUM)

        for samRun in samRuns:
            # first round of processing the sample
            if not self.getProperty("Sum").value and samRun > 0:
                self._info = None
                returned = self._focusChunks(samRun, SUFFIX, timeFilterWall, calib, splitwksp=self._splitws,\
                                             normalisebycurrent=self._normalisebycurrent,
                                             preserveEvents=preserveEvents)

                if isinstance(returned, list):
                    # Returned with a list of workspaces
                    focusedwksplist = returned
                    irun = 0
                    for run in focusedwksplist:
                        if run is not None:
                            samwksplist.append(run)
                            workspacelist.append(str(run))
                        else:
                            self.log().warning("Found a None entry in returned focused workspaces.  Index = %d." % (irun))
                        # ENDIF
                        irun += 1
                    # ENDFOR
                else:
                    run = returned
                    samwksplist.append(run)
                    workspacelist.append(str(run))
                # ENDIF
            # ENDIF
        # ENDFOR

        for (samRunIndex, samRun) in enumerate(samwksplist):
            samRun = mtd[str(samRun)]
            try:
                self.log().information("Sample Run %s:  starting number of events = %d" % (str(samRun), samRun.getNumberEvents()))
            except Exception as e:
                self.log().information("Unable to get number of events of sample run %s.  Error message: %s" % (str(samRun), str(e)))

            # Get run number
            self._info = self._getinfo(samRun.name())

            # process the container
            canRuns = self._info["container"].value
            if noRunSpecified(canRuns):
                canRun = None
            else:
                if self.getProperty("FilterCharacterizations").value:
                    canFilterWall = timeFilterWall
                else:
                    canFilterWall = (0., 0.)

                if len(canRuns) == 1:
                    canRun = canRuns[0]
                else:
                    canRun = canRuns[samRunIndex]
                if "%s_%d" % (self._instrument, canRun) in mtd:
                    canRun = mtd["%s_%d" % (self._instrument, canRun)]
                    canRun = api.ConvertUnits(InputWorkspace=canRun, OutputWorkspace=canRun, Target="TOF")
                else:
                    if self.getProperty("Sum").value:
                        canRun = self._focusAndSum(canRuns, SUFFIX, canFilterWall, calib,\
                               preserveEvents=preserveEvents)
                    else:
                        canRun = self._focusChunks(canRun, SUFFIX, canFilterWall, calib,\
                                                   normalisebycurrent=self._normalisebycurrent,
                                                   preserveEvents=preserveEvents)
                    canRun = api.ConvertUnits(InputWorkspace=canRun, OutputWorkspace=canRun, Target="TOF")
                    smoothParams = self.getProperty("BackgroundSmoothParams").value
                    if smoothParams != None and len(smoothParams)>0:
                        canRun = api.FFTSmooth(InputWorkspace=canRun, OutputWorkspace=canRun, Filter="Butterworth",\
                                               Params=smoothParams,IgnoreXBins=True,AllSpectra=True)
                workspacelist.append(str(canRun))

            # process the vanadium run
            van_run_number_list = self._info["vanadium"].value
            van_specified = not noRunSpecified(van_run_number_list)
            if van_specified:
                # get the right van run number to this sample
                if len(van_run_number_list) == 1:
                    van_run_number = van_run_number_list[0]
                else:
                    van_run_number = van_run_number_list[samRunIndex]
                # set up filter wall for van run
                if self.getProperty("FilterCharacterizations").value:
                    vanFilterWall = {'FilterByTimeStart': timeFilterWall[0], 'FilterByTimeStop': timeFilterWall[1]}
                else:
                    vanFilterWall = {'FilterByTimeStart': Property.EMPTY_DBL, 'FilterByTimeStop': Property.EMPTY_DBL}
                # get handle on workspace of this van run and make sure its unit is T.O.F
                if "%s_%d" % (self._instrument, van_run_number) in mtd:
                    # use the existing vanadium
                    van_run_ws_name = "%s_%d" % (self._instrument, van_run_number)
                    van_run_ws = mtd[van_run_ws_name]
                    van_run_ws = api.ConvertUnits(InputWorkspace=van_run_ws_name,
                                                  OutputWorkspace=van_run_ws_name,
                                                  Target="TOF")
                else:
                    # load the vanadium
                    van_run_ws_name = "%s_%d" % (self._instrument, van_run_number)
                    if self.getProperty("Sum").value:
                        van_run_ws = self._loadAndSum(van_run_number_list, van_run_ws_name, **vanFilterWall)
                    else:
                        van_run_ws = self._loadAndSum([van_run_number], van_run_ws_name, **vanFilterWall)
                    assert van_run_ws is not None

                    # load the vanadium background (if appropriate)
                    van_bkgd_run_number_list = self._info["empty"].value
                    if not noRunSpecified(van_bkgd_run_number_list):
                        if len(van_bkgd_run_number_list) == 1:
                            van_bkgd_run_number = van_bkgd_run_number_list[0]
                        else:
                            van_bkgd_run_number = van_bkgd_run_number_list[samRunIndex]
                        van_bkgd_ws_name = "%s_%d" % (self._instrument, van_bkgd_run_number)
                        if self.getProperty("Sum").value:
                            van_bkgd_ws = self._loadAndSum(van_bkgd_run_number_list, van_bkgd_ws_name, **vanFilterWall)
                        else:
                            van_bkgd_ws = self._loadAndSum([van_bkgd_run_number], van_bkgd_ws_name, **vanFilterWall)

                        if van_bkgd_ws.id() == EVENT_WORKSPACE_ID and van_bkgd_ws.getNumberEvents() <= 0:
                            # skip if background run is empty
                            pass
                        else:
                            van_run_ws = api.Minus(LHSWorkspace=van_run_ws_name,
                                                   RHSWorkspace=van_bkgd_ws_name,
                                                   OutputWorkspace=van_run_ws_name,
                                                   ClearRHSWorkspace=allEventWorkspaces(van_run_ws_name, van_bkgd_ws_name))
                        api.DeleteWorkspace(Workspace=van_bkgd_ws_name)

                    # compress events
                    if van_run_ws.id() == EVENT_WORKSPACE_ID:
                        van_run_number = api.CompressEvents(InputWorkspace=van_run_ws_name,
                                                            OutputWorkspace=van_run_ws_name,
                                                            Tolerance=self.COMPRESS_TOL_TOF)  # 10ns

                    # do the absorption correction
                    van_run_number = api.ConvertUnits(InputWorkspace=van_run_ws_name,
                                                      OutputWorkspace=van_run_ws_name,
                                                      Target="Wavelength")
                    api.SetSampleMaterial(InputWorkspace=van_run_number, ChemicalFormula="V", SampleNumberDensity=0.0721)
                    van_run_number = api.MultipleScatteringCylinderAbsorption(InputWorkspace=van_run_number, OutputWorkspace=van_run_number,
                                                                      CylinderSampleRadius=self._vanRadius)
                    van_run_number = api.ConvertUnits(InputWorkspace=van_run_number, OutputWorkspace=van_run_number, Target="TOF")

                    # focus the data
                    van_run_number = api.AlignAndFocusPowder(InputWorkspace=van_run_number, OutputWorkspace=van_run_number, CalFileName=calib,
                                                     Params=self._binning, ResampleX=self._resampleX, Dspacing=self._bin_in_dspace,
                                                     RemovePromptPulseWidth=self._removePromptPulseWidth,
                                                     CompressTolerance=self.COMPRESS_TOL_TOF,
                                                     UnwrapRef=self._LRef, LowResRef=self._DIFCref,
                                                     LowResSpectrumOffset=self._lowResTOFoffset,
                                                     CropWavelengthMin=self._wavelengthMin,
                                                     CropWavelengthMax=self._wavelengthMax,
                                                     ReductionProperties="__snspowderreduction", **(focuspos))


                    # strip peaks
                    if self.getProperty("StripVanadiumPeaks").value:
                        van_run_number = api.ConvertUnits(InputWorkspace=van_run_number, OutputWorkspace=van_run_number, Target="dSpacing")
                        van_run_number = api.StripVanadiumPeaks(InputWorkspace=van_run_number, OutputWorkspace=van_run_number, FWHM=self._vanPeakFWHM,\
                                           PeakPositionTolerance=self.getProperty("VanadiumPeakTol").value,\
                                           BackgroundType="Quadratic", HighBackground=True)
                    else:
                        self.log().information("Not strip vanadium peaks")
                    van_run_number = api.ConvertUnits(InputWorkspace=van_run_number, OutputWorkspace=van_run_number, Target="TOF")
                    van_run_number = api.FFTSmooth(InputWorkspace=van_run_number, OutputWorkspace=van_run_number, Filter="Butterworth",\
                              Params=self._vanSmoothing,IgnoreXBins=True,AllSpectra=True)
                    van_run_number = api.SetUncertainties(InputWorkspace=van_run_number, OutputWorkspace=van_run_number)
                    van_run_number = api.ConvertUnits(InputWorkspace=van_run_number, OutputWorkspace=van_run_number, Target="TOF")
                workspacelist.append(str(van_run_number))
            else:
                van_run_number = None

            if mpiRank > 0:
                return
            if samRun == 0:
                return
            # the final bit of math
            if canRun is not None:
                # must convert the sample to a matrix workspace if the can run isn't one
                if canRun.id() != EVENT_WORKSPACE_ID and samRun.id() == EVENT_WORKSPACE_ID:
                    samRun = api.ConvertToMatrixWorkspace(InputWorkspace=samRun, OutputWorkspace=samRun)
                samRun = api.Minus(LHSWorkspace=samRun, RHSWorkspace=canRun, OutputWorkspace=samRun)
                if samRun.id() == EVENT_WORKSPACE_ID:
                    samRun = api.CompressEvents(InputWorkspace=samRun, OutputWorkspace=samRun,\
                               Tolerance=self.COMPRESS_TOL_TOF) # 10ns
                canRun = str(canRun)
            if van_run_number is not None:
                samRun = api.Divide(LHSWorkspace=samRun, RHSWorkspace=van_run_number, OutputWorkspace=samRun)
                normalized = True
                samRun.getRun()['van_number'] = van_run_number.getRun()['run_number'].value
                van_run_number = str(van_run_number)
            else:
                normalized = False

            if samRun.id() == EVENT_WORKSPACE_ID:
                samRun = api.CompressEvents(InputWorkspace=samRun, OutputWorkspace=samRun,\
                           Tolerance=self.COMPRESS_TOL_TOF) # 5ns/

            # make sure there are no negative values - gsas hates them
            if self.getProperty("PushDataPositive").value != "None":
                addMin = (self.getProperty("PushDataPositive").value == "AddMinimum")
                samRun = api.ResetNegatives(InputWorkspace=samRun, OutputWorkspace=samRun, AddMinimum=addMin, ResetValue=0.)

            # write out the files
            if mpiRank == 0:
                if self._scaleFactor != 1.:
                    samRun = api.Scale(samRun, Factor=self._scaleFactor, OutputWorkspace=samRun)
                self._save(samRun, self._info, normalized, False)
                samRun = str(samRun)
            #mtd.releaseFreeMemory()

        # ENDFOR

        # convert everything into d-spacing
        workspacelist = set(workspacelist) # only do each workspace once
        if HAVE_MPI is False:
            for wksp in workspacelist:
                wksp = api.ConvertUnits(InputWorkspace=wksp, OutputWorkspace=wksp, Target=self.getProperty("FinalDataUnits").value)

        return

    def _loadCharacterizations(self):
        self._focusPos = {}
        charFilename = self.getProperty("CharacterizationRunsFile").value
        expIniFilename = self.getProperty("ExpIniFilename").value

        if charFilename is None or len(charFilename) <= 0:
            self.iparmFile = None
            return

        results = api.PDLoadCharacterizations(Filename=charFilename,
                                              ExpIniFilename=expIniFilename,
                                              OutputWorkspace="characterizations")
        # export the characterizations table
        self._charTable = results[0]
        self.declareProperty(ITableWorkspaceProperty("CharacterizationsTable", "characterizations", Direction.Output))
        self.setProperty("CharacterizationsTable", self._charTable)

        # get the focus positions from the properties
        self.iparmFile = results[1]
        self._focusPos['PrimaryFlightPath'] = results[2]
        self._focusPos['SpectrumIDs'] = results[3]
        self._focusPos['L2'] = results[4]
        self._focusPos['Polar'] = results[5]
        self._focusPos['Azimuthal'] = results[6]

    #pylint: disable=too-many-branches
    def _loadData(self, run_number, extension, filterWall=None, out_ws_name=None, **chunk):
        """ Load data optionally by chunk strategy
        Purpose:
            Load a complete or partial run, filter bad pulses.
            Output workspace name is formed as
            - user specified (highest priority)
            - instrument-name_run-number_0 (no chunking)
            - instrument-name_run-number_X: X is either ChunkNumber of SpectrumMin
        Requirements:
            1. run number is integer and larger than 0
        Guarantees:
            A workspace is created with name described above
        :param run_number:
        :param extension: file extension
        :param filterWall:
        :param out_ws_name: name of output workspace specified by user. it will override the automatic name
        :param chunk:
        :return:
        """
        # Check requirements
        assert isinstance(run_number, int), 'Input run number must be integer but not %s.' % str(type(run_number))
        assert run_number > 0, 'Input run number must be larger than 0 but not %d.' % run_number
        assert (chunk is None) or isinstance(chunk, dict), 'Input chunk must be either a dictionary or None.'

        base_name = "%s_%d" % (self._instrument, run_number)
        filename = base_name + extension
        # give out the default output workspace name
        if out_ws_name is None:
            if chunk:
                if "ChunkNumber" in chunk:
                    out_ws_name = base_name + "__chk%d" % (int(chunk["ChunkNumber"]))
                elif "SpectrumMin" in chunk:
                    seq_number = 1 + int(chunk["SpectrumMin"])/(int(chunk["SpectrumMax"])-int(chunk["SpectrumMin"]))
                    out_ws_name = base_name + "__chk%d" % seq_number
            else:
                out_ws_name = "%s_%d" % (base_name, 0)
        # END-IF

        # Specify the other chunk information including Precount, FilterByTimeStart and FilterByTimeStop.
        if extension.endswith("_event.nxs") or extension.endswith(".nxs.h5"):
            chunk["Precount"] = True
            if filterWall is not None:
                if filterWall[0] > 0.:
                    chunk["FilterByTimeStart"] = filterWall[0]
                if filterWall[1] > 0.:
                    chunk["FilterByTimeStop"] = filterWall[1]

        # Call Mantid's Load algorithm to load complete or partial data
        wksp = api.Load(Filename=filename, OutputWorkspace=out_ws_name, **chunk)
        assert wksp is not None, 'Mantid Load algorithm does not return a workspace.'

        # Log output
        try:
            self.log().debug("Load run %s: number of events = %d" % (str(run_number), wksp.getNumberEvents()))
        except Exception as e:
            self.log().debug("Load run %s: unable to get events of %s.  Error message: %s" % (str(run_number),                                                                                     str(wksp), str(e)))
        if HAVE_MPI:
            msg = "MPI Task = %s ;" % (str(mpiRank))
            try:
                msg += "Number Events = " + str(wksp.getNumberEvents())
            except Exception as e:
                msg += "Unable to get events of %s.  Error message: %s" % (str(wksp), str(e))
            self.log().debug(msg)

        # filter bad pulses
        if self._filterBadPulses > 0.:
            # record number of events of the original workspace
            is_event_ws = isinstance(wksp, mantid.api.IEventWorkspace)
            if is_event_ws is True:
                # Event workspace: record original number of events
                num_original_events = wksp.getNumberEvents()
            else:
                num_original_events = -1

            # filter bad pulse
            wksp = api.FilterBadPulses(InputWorkspace=out_ws_name, OutputWorkspace=out_ws_name,
                                       LowerCutoff=self._filterBadPulses)
            assert wksp is not None, 'Returned value from FilterBadPulses cannot be None'

            if is_event_ws is True:
                # Event workspace
                message = "FilterBadPulses reduces number of events from %d to %d (under %s percent) " \
                          "of workspace %s." % (num_original_events, wksp.getNumberEvents(),
                                                str(self._filterBadPulses), str(wksp))
                self.log().information(message)
        # END-IF (filter bad pulse)

        return wksp

    def _getStrategy(self, runnumber, extension):
        """
        Get chunking strategy by calling mantid algorithm 'DetermineChunking'
        :param runnumber:
        :param extension:
        :return: a list of dictionary.  Each dictionary is a row in table workspace
        """
        # generate the workspace name
        assert isinstance(runnumber, int)
        file_name = "%s_%d" % (self._instrument, runnumber) + extension
        assert os.path.exists(file_name), 'NeXus file %s does not exist.' % file_name

        self.log().debug("[Fx116] Run file Name : %s,\t\tMax chunk size: %s" % (file_name, str(self._chunks)))
        chunks = api.DetermineChunking(Filename=file_name, MaxChunkSize=self._chunks)

        strategy = []
        for row in chunks:
            strategy.append(row)

        # For table with no rows
        if len(strategy) == 0:
            strategy.append({})

        # delete chunks workspace
        chunks = str(chunks)
        mtd.remove(chunks)

        return strategy

    def __logChunkInfo(self, chunk):
        keys = chunk.keys()
        keys.sort()

        keys = [ str(key) + "=" + str(chunk[key]) for key in keys ]
        self.log().information("Working on chunk [" + ", ".join(keys) + "]")

    def checkInfoMatch(self, left, right):
        if (left["frequency"].value is not None) and (right["frequency"].value is not None) \
           and (abs(left["frequency"].value - right["frequency"].value)/left["frequency"].value > .05):
            raise RuntimeError("Cannot add incompatible frequencies (%f!=%f)" \
                               % (left["frequency"].value, right["frequency"].value))
        if (left["wavelength"].value is not None) and (right["wavelength"].value is not None) \
                   and abs(left["wavelength"].value - right["wavelength"].value)/left["wavelength"].value > .05:
            raise RuntimeError("Cannot add incompatible wavelengths (%f != %f)" \
                               % (left["wavelength"].value, right["wavelength"].value))

    def _loadAndSum(self, runnumbers, outName, **filterWall):
        names=["%s_%d" % (self._instrument, runNum) for runNum in runnumbers]

        sumRun = None
        info = None
        SUFFIX = self.getProperty("Extension").value

        for name in names:
            self.log().information("[Sum] processing %s" % name)
            temp = api.LoadEventAndCompress(Filename=name+SUFFIX, OutputWorkspace=name,
                                            MaxChunkSize=self._chunks, FilterBadPulses=self._filterBadPulses,
                                            CompressTOFTolerance=self.COMPRESS_TOL_TOF,
                                            **filterWall)
            tempinfo = self._getinfo(temp)
            if sumRun is None:
                sumRun = temp
                info = tempinfo
            else:
                self.checkInfoMatch(info, tempinfo)

                sumRun = api.Plus(LHSWorkspace=sumRun, RHSWorkspace=temp, OutputWorkspace=sumRun,
                                  ClearRHSWorkspace=allEventWorkspaces(sumRun, temp))
                if sumRun.id() == EVENT_WORKSPACE_ID:
                    sumRun = api.CompressEvents(InputWorkspace=sumRun, OutputWorkspace=sumRun,\
                                                Tolerance=self.COMPRESS_TOL_TOF) # 10ns
                api.DeleteWorkspace(str(temp))

        if str(sumRun) != outName:
            sumRun = api.RenameWorkspace(InputWorkspace=sumRun, OutputWorkspace=outName)

        try:
            if self._normalisebycurrent is True:
                sumRun = api.NormaliseByCurrent(InputWorkspace=temp,
                                              OutputWorkspace=temp)
                temp.getRun()['gsas_monitor'] = 1
        except Exception, e:
            self.log().warning(str(e))

        return sumRun


    #pylint: disable=too-many-arguments
    def _focusAndSum(self, runnumbers, extension, filterWall, calib, preserveEvents=True):
        """Load, sum, and focus data in chunks"""
        sumRun = None
        info = None

        for temp in runnumbers:
            runnumber = temp
            self.log().information("[Sum] Process run number %s. " %(str(runnumber)))

            temp = self._focusChunks(temp, extension, filterWall, calib,\
                                     normalisebycurrent=False,
                                     preserveEvents=preserveEvents)
            tempinfo = self._getinfo(temp)

            if sumRun is None:
                sumRun = temp
                info = tempinfo
            else:
                self.checkInfoMatch(info, tempinfo)

                sumRun = api.Plus(LHSWorkspace=sumRun, RHSWorkspace=temp, OutputWorkspace=sumRun,
                                  ClearRHSWorkspace=allEventWorkspaces(sumRun, temp))
                if sumRun.id() == EVENT_WORKSPACE_ID:
                    sumRun = api.CompressEvents(InputWorkspace=sumRun, OutputWorkspace=sumRun,\
                                                Tolerance=self.COMPRESS_TOL_TOF) # 10ns
                api.DeleteWorkspace(str(temp))
            # ENDIF
        # ENDFOR (processing each)
        if self._normalisebycurrent is True:
            sumRun = api.NormaliseByCurrent(InputWorkspace=sumRun,
                                            OutputWorkspace=sumRun)
            sumRun.getRun()['gsas_monitor'] = 1

        return sumRun


    #pylint: disable=too-many-arguments,too-many-locals,too-many-branches
    def _focusChunks(self, runnumber, extension, filterWall, calib,
                     normalisebycurrent, splitwksp=None, preserveEvents=True):
        """ Load, (optional) split and focus data in chunks

        Arguments:
         - runnumber : integer for run number
         - normalisebycurrent: Set to False if summing runs for correct math
         - splitwksp:  SplittersWorkspace (if None then no split)
         - filterWall: Enabled if splitwksp is defined

        Return:
        """
        # generate the workspace name
        wksp = "%s_%d" % (self._instrument, runnumber)
        self.log().information("_focusChunks(): runnumber = %d, extension = %s" % (runnumber, extension))

        # get chunk strategy for parallel processing (MPI)
        strategy = self._getStrategy(runnumber, extension)

        # determine event splitting by checking filterWall and number of output workspaces from _focusChunk
        do_split_raw_wksp, num_out_wksp = self._determine_workspace_splitting(splitwksp, filterWall)

        # Set up the data structure to hold and control output workspaces
        output_wksp_list = [None] * num_out_wksp
        is_first_chunk_list = [True] * num_out_wksp

        self.log().debug("F1141A: Number of workspace to process = %d" % num_out_wksp)

        # reduce data by chunks
        chunk_index = -1
        # FIXME/TODO/Now: there are 2 places that determine the base name
        base_name = "%s_%d" % (self._instrument, runnumber)
        for chunk in strategy:
            # progress on chunk index
            chunk_index += 1
            # Load chunk, i.e., partial data into Mantid
            raw_ws_chunk = self._loadData(runnumber, extension, filterWall, **chunk)
            raw_ws_name_chunk = raw_ws_chunk.name()

            if self._info is None:
                self._info = self._getinfo(raw_ws_name_chunk)

            # Log information for current chunk
            self.__logChunkInfo(chunk)
            if raw_ws_chunk.id() == EVENT_WORKSPACE_ID:
                # Event workspace
                self.log().debug("F1141C: There are %d events after data is loaded in workspace %s." % (
                    raw_ws_chunk.getNumberEvents(), raw_ws_name_chunk))

            # Split the workspace if it is required
            if do_split_raw_wksp is True:
                output_ws_name_list = self._split_workspace(raw_ws_name_chunk, splitwksp)
            else:
                # Histogram data
                output_ws_name_list = [raw_ws_name_chunk]
            # ENDIF
            # check
            if num_out_wksp != len(output_ws_name_list):
                self.log().warning('Projected number of output workspaces %d must be same as '
                                   'that of real output workspaces %d.' % (num_out_wksp, len(output_ws_name_list)))
                num_out_wksp = output_ws_name_list

            # log
            msg = "[Fx1142] Workspace of chunk %d is %d (vs. estimated %d). \n" % (
                chunk_index, len(output_ws_name_list), num_out_wksp)
            for iws in xrange(len(output_ws_name_list)):
                ws = output_ws_name_list[iws]
                msg += "%s\t\t" % (str(ws))
                if iws %5 == 4:
                    msg += "\n"
            self.log().debug(msg)

            # Do align and focus
            num_out_wksp = len(output_ws_name_list)
            for split_index in xrange(num_out_wksp):
                # Get workspace name
                out_ws_name_chunk_split = output_ws_name_list[split_index]
                # Align and focus
                # focuspos = self._focusPos
                self.log().notice('Align and focus workspace %s' % out_ws_name_chunk_split)
                out_ws_c_s = api.AlignAndFocusPowder(InputWorkspace=out_ws_name_chunk_split,
                                                     OutputWorkspace=out_ws_name_chunk_split,
                                                     CalFileName=calib,
                                                     Params=self._binning,
                                                     ResampleX=self._resampleX,
                                                     Dspacing=self._bin_in_dspace,
                                                     PreserveEvents=preserveEvents,
                                                     RemovePromptPulseWidth=self._removePromptPulseWidth,
                                                     CompressTolerance=self.COMPRESS_TOL_TOF,
                                                     UnwrapRef=self._LRef,
                                                     LowResRef=self._DIFCref,
                                                     LowResSpectrumOffset=self._lowResTOFoffset,
                                                     CropWavelengthMin=self._wavelengthMin,
                                                     CropWavelengthMax=self._wavelengthMax,
                                                     ReductionProperties="__snspowderreduction",
                                                     **self._focusPos)
                assert out_ws_c_s is not None
                # logging (ignorable)
                for iws in xrange(out_ws_c_s.getNumberHistograms()):
                    spec = out_ws_c_s.getSpectrum(iws)
                    self.log().debug("[DBx131] ws %d: spectrum ID = %d. " % (iws, spec.getSpectrumNo()))
                if out_ws_c_s.id() == EVENT_WORKSPACE_ID:
                    self.log().information('After being aligned and focused, workspace %s: Number of events = %d '
                                           'of chunk %d ' % (out_ws_c_s.name(), out_ws_c_s.getNumberEvents(),
                                                             chunk_index))
            # END-FOR-Splits

            # Merge among chunks
            for split_index in range(num_out_wksp):

                # determine the final workspace name
                final_out_ws_name = base_name
                if num_out_wksp > 1:
                    final_out_ws_name += '_%d' % split_index

                if is_first_chunk_list[split_index] is True:
                    # Rename if it is the first chunk that is finished
                    self.log().debug("[F1145] Slot %d is renamed to %s" % (split_index, final_out_ws_name))
                    temp_out_ws = api.RenameWorkspace(InputWorkspace=output_ws_name_list[split_index],
                                                      OutputWorkspace=final_out_ws_name)
                    assert temp_out_ws is not None

                    output_wksp_list[split_index] = final_out_ws_name
                    is_first_chunk_list[split_index] = False
                else:
                    # Add this chunk of workspace to the final one
                    clear_rhs_ws = allEventWorkspaces(output_wksp_list[split_index], output_ws_name_list[split_index])
                    temp_out_ws = api.Plus(LHSWorkspace=output_wksp_list[split_index],
                                           RHSWorkspace=output_ws_name_list[split_index],
                                           OutputWorkspace=output_wksp_list[split_index],
                                           ClearRHSWorkspace=clear_rhs_ws)
                    assert temp_out_ws is not None

                    # Delete the chunk workspace
                    api.DeleteWorkspace(output_ws_name_list[split_index])
                # END-IF-ELSE
        # END-FOR-Chunks

        # Sum workspaces for all mpi tasks
        if HAVE_MPI:
            for split_index in xrange(num_out_wksp):
                output_wksp_list[split_index] = api.GatherWorkspaces(
                    InputWorkspace=output_wksp_list[split_index],
                    PreserveEvents=preserveEvents,
                    AccumulationMethod="Add",
                    OutputWorkspace=output_wksp_list[split_index])
        # ENDIF MPI

        if self._chunks > 0:
            # When chunks are added, proton charge is summed for all chunks
            for split_index in xrange(num_out_wksp):
                output_wksp_list[split_index].getRun().integrateProtonCharge()
        # ENDIF

        if (self.iparmFile is not None) and (len(self.iparmFile) > 0):
            # When chunks are added, add iparamFile
            for split_index in xrange(num_out_wksp):
                output_wksp_list[split_index].getRun()['iparm_file'] = self.iparmFile

        # Compress events
        for split_index in xrange(num_out_wksp):
            temp_out_ws = self.get_workspace(output_wksp_list[split_index])
            if temp_out_ws.id() == EVENT_WORKSPACE_ID:
                temp_out_ws = api.CompressEvents(InputWorkspace=output_wksp_list[split_index],
                                                 OutputWorkspace=output_wksp_list[split_index],
                                                 Tolerance=self.COMPRESS_TOL_TOF) # 100ns
                assert temp_out_ws is not None

            try:
                if normalisebycurrent is True:
                    output_wksp_list[split_index] = api.NormaliseByCurrent(InputWorkspace=output_wksp_list[split_index],
                                                                           OutputWorkspace=output_wksp_list[split_index])
                    output_wksp_list[split_index].getRun()['gsas_monitor'] = 1
            except Exception, e:
                self.log().warning(str(e))

            propertyName = "OutputWorkspace%s" % str(output_wksp_list[split_index])
            self.log().warning(propertyName)
            self.declareProperty(WorkspaceProperty(propertyName, str(output_wksp_list[split_index]), Direction.Output))
            self.setProperty(propertyName, output_wksp_list[split_index])
            self._save(output_wksp_list[split_index], self._info, False, True)
            self.log().information("Done focussing data of %d." % (split_index))

        self.log().information("[E1207] Number of workspace in workspace list after clean = %d. " %(len(output_wksp_list)))

        # About return
        if splitwksp is None:
            return output_wksp_list[0]
        else:
            return output_wksp_list

    def _getinfo(self, wksp_name):
        """ Get characterization information of a certain workspace
        Purpose:

        Requirements:
        1. wksp_name is string and AnalysisDataService has such workspace
        Guarantees:

        :param wksp_name:
        :return:
        """
        # Check requirements
        assert isinstance(wksp_name, str)
        assert self.does_workspace_exist(wksp_name)

        # Determine characterization
        if mtd.doesExist("characterizations"):
            # get the correct row of the table if table workspace 'charactersizations' exists

            #pylint: disable=unused-variable
            charac = api.PDDetermineCharacterizations(InputWorkspace=wksp_name,
                                                      Characterizations="characterizations",
                                                      ReductionProperties="__snspowderreduction",
                                                      BackRun=self.getProperty("BackgroundNumber").value,
                                                      NormRun=self.getProperty("VanadiumNumber").value,
                                                      NormBackRun=self.getProperty("VanadiumBackgroundNumber").value,
                                                      FrequencyLogNames=self.getProperty("FrequencyLogNames").value,
                                                      WaveLengthLogNames=self.getProperty("WaveLengthLogNames").value)
        else:
            charac = api.PDDetermineCharacterizations(InputWorkspace=wksp_name,
                                                      ReductionProperties="__snspowderreduction",
                                                      BackRun=self.getProperty("BackgroundNumber").value,
                                                      NormRun=self.getProperty("VanadiumNumber").value,
                                                      NormBackRun=self.getProperty("VanadiumBackgroundNumber").value,
                                                      FrequencyLogNames=self.getProperty("FrequencyLogNames").value,
                                                      WaveLengthLogNames=self.getProperty("WaveLengthLogNames").value)

        # convert the result into a dict
        return PropertyManagerDataService.retrieve("__snspowderreduction")

    def _save(self, wksp, info, normalized, pdfgetn):
        prefix = str(wksp)
        if len(self._outPrefix) > 0: # non-empty string
            prefix = self._outPrefix
        filename = os.path.join(self._outDir, prefix)
        if pdfgetn:
            self.log().notice("Saving 'pdfgetn' is deprecated. Use PDtoPDFgetN instead.")
            if "pdfgetn" in self._outTypes:
                pdfwksp = str(wksp)+"_norm"
                pdfwksp = api.SetUncertainties(InputWorkspace=wksp, OutputWorkspace=pdfwksp, SetError="sqrt")
                api.SaveGSS(InputWorkspace=pdfwksp, Filename=filename+".getn", SplitFiles=False, Append=False,\
                        MultiplyByBinWidth=False, Bank=info["bank"].value, Format="SLOG", ExtendedHeader=True)
                api.DeleteWorkspace(pdfwksp)
            return # don't do the other bits of saving
        if "gsas" in self._outTypes:
            api.SaveGSS(InputWorkspace=wksp, Filename=filename+".gsa", SplitFiles=False, Append=False,\
                    MultiplyByBinWidth=normalized, Bank=info["bank"].value, Format="SLOG", ExtendedHeader=True)
        if "fullprof" in self._outTypes:
            api.SaveFocusedXYE(InputWorkspace=wksp, StartAtBankNumber=info["bank"].value, Filename=filename+".dat")
        if "topas" in self._outTypes:
            api.SaveFocusedXYE(InputWorkspace=wksp, StartAtBankNumber=info["bank"].value, Filename=filename+".xye",
                               Format="TOPAS")
        if "nexus" in self._outTypes:
            api.ConvertUnits(InputWorkspace=wksp, OutputWorkspace=wksp, Target=self.getProperty("FinalDataUnits").value)
            #api.Rebin(InputWorkspace=wksp, OutputWorkspace=wksp, Params=self._binning) # crop edges
            api.SaveNexus(InputWorkspace=wksp, Filename=filename+".nxs")

        # always save python script - this is broken because the history isn't
        # attached until the algorithm is finished
        api.GeneratePythonScript(InputWorkspace=wksp, Filename=filename+".py")

        return

    def _getTimeFilterWall(self, splitws, samrun, extension):
        """ Get filter wall from splitter workspace, i.e.,
        get the earlies and latest TIME stamp in input splitter workspace

        Arguments:
         - splitws      : splitters workspace
         - runstarttime : total nanoseconds of run start time (Mantid DateAndTime)

        Return: tuple of start-time and stop-time relative to run start time and in unit of second
                If there is no split workspace defined, filter is (0., 0.) as the default
        """
        # None case
        if splitws is None:
            self.log().warning("Split workspace is None.  Unable to make a filter wall.  Return with default value. ")
            return (0.0, 0.0)

        # Load data
        name = "%s_%d" % (self._instrument, samrun)
        filename = name + extension
        metawsname = "temp_"+name

        metawksp = api.Load(Filename=str(filename), OutputWorkspace=str(metawsname), MetaDataOnly=True)
        if metawksp is None:
            self.log().warning("Unable to open file %s" % (filename))
            return (0.0, 0.0)

        # Get start time
        runstarttimens = metawksp.getRun().startTime().totalNanoseconds()

        numrow = splitws.rowCount()

        # Searching for the
        tmin_absns = splitws.cell(0,0)
        tmax_absns = splitws.cell(0,1)

        for r in xrange(1, numrow):
            timestart = splitws.cell(r, 0)
            timeend = splitws.cell(r, 1)
            if timestart < tmin_absns:
                tmin_absns = timestart
            if timeend > tmax_absns:
                tmax_absns = timeend
        # ENDFOR

        tmin = (tmin_absns - runstarttimens) * 1.0E-9
        tmax = (tmax_absns - runstarttimens) * 1.0E-9

        filterWall = (tmin, tmax)

        api.DeleteWorkspace(Workspace=metawsname)

        return filterWall

    def getNumberOfSplittedWorkspace(self, splitwksp):
        """ Get number of splitted workspaces due to input splitwksp

        Return : integer
        """
        # splitws = mtd["PG3_9829_event_splitters"]
        splitws = AnalysisDataService.retrieve(str(splitwksp))
        numrows = splitws.rowCount()
        wscountdict = {}
        for r in xrange(numrows):
            wsindex = splitws.cell(r,2)
            wscountdict[wsindex] = 0

        return len(wscountdict.keys())

    @staticmethod
    def does_workspace_exist(workspace_name):
        """
        Purpose: Check whether a workspace exists in AnalysisDataService
        :param workspace_name:
        :return:
        """
        assert isinstance(workspace_name, str)

        return AnalysisDataService.doesExist(workspace_name)

    def get_workspace(self, workspace_name):
        """
        Purpose: Get the reference of a workspace
        Requirements:
            1. workspace_name is a string
        Guarantees:
            The reference of the workspace with same is returned.
        :exception RuntimeError: a RuntimeError from Mantid will be thrown
        :param workspace_name:
        :return:
        """
        assert isinstance(workspace_name, str)

        return AnalysisDataService.retrieve(workspace_name)

    def _determine_workspace_splitting(self, split_wksp, filter_wall):
        """
        :return: do_split_raw_wksp, num_out_wksp
        """
        do_split_raw_wksp = False
        num_out_wksp = 1

        if split_wksp is not None:
            # Use splitting workspace

            # Check consistency with filterWall
            if filter_wall[0] < 1.0E-20 and filter_wall[1] < 1.0E-20:
                # Default definition of filterWall when there is no split workspace specified.
                raise RuntimeError("It is impossible to have a not-NONE splitters workspace and (0,0) time filter wall.")
            # ENDIF

            # Note: Unfiltered workspace (remainder) is not considered here
            num_out_wksp = self.getNumberOfSplittedWorkspace(split_wksp)
            num_splitters = split_wksp.rowCount()

            # Do explicit FilterEvents if number of splitters is larger than 1.
            # If number of splitters is equal to 1, then filterWall will do the job itself.
            if num_splitters > 1:
                do_split_raw_wksp = True
            self.log().debug("[Fx948] Number of split workspaces = %d; Do split = %s" % (num_out_wksp, str(do_split_raw_wksp)))
        # ENDIF

        return do_split_raw_wksp, num_out_wksp

    def _split_workspace(self, raw_ws_name, split_ws_name):
        """ Split workspace
        Purpose:
            Split a workspace
        Requirements:
            1. raw_ws_name is a string
            2. an event workspace with this name exists
        Guarantees:
            Raw input workspace is split
        :param raw_ws_name:
        :param split_ws_name:
        :return: list of strings as output workspaces
        """
        # Check requirements
        assert isinstance(raw_ws_name, str), 'Raw workspace name must be a string.'
        assert isinstance(split_ws_name, str)
        assert self.does_workspace_exist(split_ws_name)

        raw_ws = self.get_workspace(workspace_name=raw_ws_name)
        assert raw_ws.id() == EVENT_WORKSPACE_ID, 'Input workspace for splitting must be an EventWorkspace.'

        # Splitting workspace
        self.log().information("SplitterWorkspace = %s, Information Workspace = %s. " % (
                split_ws_name, str(self._splitinfotablews)))

        base_name = raw_ws_name
        if self._splitinfotablews is None:
            # split without information table
            api.FilterEvents(InputWorkspace=raw_ws_name, OutputWorkspaceBaseName=base_name,
                             SplitterWorkspace=split_ws_name, GroupWorkspaces=True)
        else:
            # split with information table
            api.FilterEvents(InputWorkspace=raw_ws_name, OutputWorkspaceBaseName=base_name,
                             SplitterWorkspace=split_ws_name, InformationWorkspace = str(self._splitinfotablews),
                             GroupWorkspaces=True)
        # ENDIF

        # Get workspace group for names of split workspace
        wsgroup = mtd[base_name]
        tempwsnamelist = wsgroup.getNames()
        # logging
        dbstr = "[Fx951] Splitted workspace names: "
        for ws_name in tempwsnamelist:
            dbstr += "%s, " % (ws_name)
        self.log().debug(dbstr)

        # Build the list of workspaces' names for return
        out_ws_name_list = []
        # FIXME Keep in mind to use this option.
        # keepremainder = self.getProperty("KeepRemainder").value
        for ws_name in tempwsnamelist:
            this_ws = self.get_workspace(ws_name)
            if ws_name.endswith("_unfiltered") is False:
                # workspace to save
                out_ws_name_list.append(ws_name)
            else:
                # events that are not excluded by filters. delete the workspace
                api.DeleteWorkspace(Workspace=this_ws)
        # END-FOR

        return out_ws_name_list

# Register algorithm with Mantid.
AlgorithmFactory.subscribe(SNSPowderReduction)
