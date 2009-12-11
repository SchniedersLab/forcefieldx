@REM ----------------------------------------------------------------------------
@REM Copyright 2001-2004 The Apache Software Foundation.
@REM
@REM Licensed under the Apache License, Version 2.0 (the "License");
@REM you may not use this file except in compliance with the License.
@REM You may obtain a copy of the License at
@REM
@REM      http://www.apache.org/licenses/LICENSE-2.0
@REM
@REM Unless required by applicable law or agreed to in writing, software
@REM distributed under the License is distributed on an "AS IS" BASIS,
@REM WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
@REM See the License for the specific language governing permissions and
@REM limitations under the License.
@REM ----------------------------------------------------------------------------
@REM

@echo off

set ERROR_CODE=0

:init
@REM Decide how to startup depending on the version of windows

@REM -- Win98ME
if NOT "%OS%"=="Windows_NT" goto Win9xArg

@REM set local scope for the variables with windows NT shell
if "%OS%"=="Windows_NT" @setlocal

@REM -- 4NT shell
if "%eval[2+2]" == "4" goto 4NTArgs

@REM -- Regular WinNT shell
set CMD_LINE_ARGS=%*
goto WinNTGetScriptDir

@REM The 4NT Shell from jp software
:4NTArgs
set CMD_LINE_ARGS=%$
goto WinNTGetScriptDir

:Win9xArg
@REM Slurp the command line arguments.  This loop allows for an unlimited number
@REM of agruments (up to the command line limit, anyway).
set CMD_LINE_ARGS=
:Win9xApp
if %1a==a goto Win9xGetScriptDir
set CMD_LINE_ARGS=%CMD_LINE_ARGS% %1
shift
goto Win9xApp

:Win9xGetScriptDir
set SAVEDIR=%CD%
%0\
cd %0\..\.. 
set BASEDIR=%CD%
cd %SAVEDIR%
set SAVE_DIR=
goto repoSetup

:WinNTGetScriptDir
set BASEDIR=%~dp0\..

:repoSetup


if "%JAVACMD%"=="" set JAVACMD=java

if "%REPO%"=="" set REPO=%BASEDIR%\..\lib

set CLASSPATH="%BASEDIR%"\etc;"%REPO%"\macosx\AppleJavaExtensions\1.3\AppleJavaExtensions-1.3.jar;"%REPO%"\commons-io\commons-io\1.4\commons-io-1.4.jar;"%REPO%"\commons-lang\commons-lang\2.4\commons-lang-2.4.jar;"%REPO%"\org\codehaus\groovy\groovy-all\1.6.5\groovy-all-1.6.5.jar;"%REPO%"\junit\junit\4.6\junit-4.6.jar;"%REPO%"\org\apache\ant\ant\1.7.1\ant-1.7.1.jar;"%REPO%"\org\apache\ant\ant-launcher\1.7.1\ant-launcher-1.7.1.jar;"%REPO%"\jline\jline\0.9.94\jline-0.9.94.jar;"%REPO%"\net\java\dev\jogl\jogl\1.1.1\jogl-1.1.1.jar;"%REPO%"\net\java\dev\jogl\gluegen-rt\1.1.1\gluegen-rt-1.1.1.jar;"%REPO%"\java3d\j3dcore\1.5.2\j3dcore-1.5.2.jar;"%REPO%"\java3d\j3dutils\1.5.2\j3dutils-1.5.2.jar;"%REPO%"\java3d\j3dvrml\1.5.2\j3dvrml-1.5.2.jar;"%REPO%"\java3d\vecmath\1.5.2\vecmath-1.5.2.jar;"%REPO%"\edu\rit\pj\pj\1.0\pj-1.0.jar;"%REPO%"\javax\help\javahelp\2.0.02\javahelp-2.0.02.jar;"%REPO%"/com/kenai/ffx/algorithms/1.0/algorithms-1.0.jar;"%REPO%"/com/kenai/ffx/crystal/1.0/crystal-1.0.jar;"%REPO%"/com/kenai/ffx/numerics/1.0/numerics-1.0.jar;"%REPO%"/com/kenai/ffx/parsers/1.0/parsers-1.0.jar;"%REPO%"/com/kenai/ffx/potentials/1.0/potentials-1.0.jar;"%REPO%"/com/kenai/ffx/ui/1.0/ui-1.0.jar

set EXTRA_JVM_ARGUMENTS=-Xms128M -Xmx256M

goto endInit

@REM Reaching here means variables are defined and arguments have been captured
:endInit

%JAVACMD% %JAVA_OPTS% %EXTRA_JVM_ARGUMENTS% -classpath %CLASSPATH_PREFIX%;%CLASSPATH% -Djava.library.path="%BASEDIR%/windows/jogl-32" -Dapp.name="ffx" -Dj3d.rend="jogl" -Dapp.repo="%REPO%" -Dbasedir="%BASEDIR%" -Djava.awt.headless="true" ffx.Main %CMD_LINE_ARGS%
if ERRORLEVEL 1 goto error
goto end

:error
if "%OS%"=="Windows_NT" @endlocal
set ERROR_CODE=1

:end
@REM set local scope for the variables with windows NT shell
if "%OS%"=="Windows_NT" goto endNT

@REM For old DOS remove the set variables from ENV - we assume they were not set
@REM before we started - at least we don't leave any baggage around
set CMD_LINE_ARGS=
goto postExec

:endNT
@endlocal

:postExec

if "%FORCE_EXIT_ON_ERROR%" == "on" (
  if %ERROR_CODE% NEQ 0 exit %ERROR_CODE%
)

exit /B %ERROR_CODE%
