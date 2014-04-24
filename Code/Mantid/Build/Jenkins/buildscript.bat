:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: WINDOWS SCRIPT TO DRIVE THE JENKINS BUILDS OF MANTID.
::
:: Notes:
::
:: WORKSPACE & JOB_NAME are environment variables that are set by Jenkins.
:: BUILD_THREADS & PARAVIEW_DIR should be set in the configuration of each slave.
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: If pvnext in JOB_NAME, use PARAVIEW_NEXT_DIR for LOCAL_PARAVIEW_DIR else 
:: use PARAVIEW_DIR
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if "%JOB_NAME%"=="%JOB_NAME:pvnext=%" (
    set LOCAL_PARAVIEW_DIR=%PARAVIEW_NEXT_DIR%
) else (
    set LOCAL_PARAVIEW_DIR=%PARAVIEW_DIR%
)

echo "A: %LOCAL_PARAVIEW_DIR%"

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Get or update the third party dependencies
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd %WORKSPACE%\Code
call fetch_Third_Party win64
cd %WORKSPACE%

set PATH=%WORKSPACE%\Code\Third_Party\lib\win64;%WORKSPACE%\Code\Third_Party\lib\win64\Python27;%LOCAL_PARAVIEW_DIR%\bin\Release;%PATH%
echo "B: %PATH%"
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Check whether this is a clean build (must have 'clean' in the job name)
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set DOC_IMAGES=
if "%JOB_NAME%"=="%JOB_NAME:clean=%" (
    set CLEANBUILD=no
) else  (
    set CLEANBUILD=yes
    rmdir /S /Q build
    if NOT "%JOB_NAME%"=="%JOB_NAME:master=%" (
        set DOC_IMAGES=-DQT_ASSISTANT_FETCH_IMAGES=ON
    )
)

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Create the build directory if it doesn't exist
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
md %WORKSPACE%\build
cd %WORKSPACE%\build

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: CMake configuration
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
"C:\Program Files (x86)\CMake 2.8\bin\cmake.exe" -G "Visual Studio 11 Win64" -DENABLE_CPACK=ON -DMAKE_VATES=ON -DParaView_DIR=%LOCAL_PARAVIEW_DIR% -DUSE_PRECOMPILED_HEADERS=ON %DOC_IMAGES% ..\Code\Mantid

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Build step
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
msbuild /nologo /m:%BUILD_THREADS% /nr:false /p:Configuration=Release Mantid.sln
if ERRORLEVEL 1 exit /B %ERRORLEVEL%

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Run the tests
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
"C:\Program Files (x86)\CMake 2.8\bin\ctest.exe" -C Release -j%BUILD_THREADS% --schedule-random --output-on-failure -E MantidPlot
:: Run GUI tests serially
ctest -C Release --output-on-failure -R MantidPlot

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Create the install kit if this is a clean build
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
if "%CLEANBUILD%" EQU "yes" (
    msbuild /nologo /m:%BUILD_THREADS% /nr:false /p:Configuration=Release docs/qtassistant/qtassistant.vcxproj
    cpack -C Release --config CPackConfig.cmake
)
