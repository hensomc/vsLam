@rem  BAT file generated from the c:/msc.software/msc_nastran_student_edition/2014/nastran/msc20140/win64i4/nastran.lcl
@rem  template file for Windows NT/2K/XP/.NET 2003/...
@rem  for Local Job processing
@rem
@rem Run MSC Nastran on the local node.
@rem
@rem THIS PROGRAM IS CONFIDENTIAL AND A TRADE SECRET OF THE MSC.SOFTWARE
@rem CORPORATION.  THE RECEIPT OR POSSESSION OF THIS PROGRAM DOES NOT CONVEY ANY
@rem RIGHTS TO REPRODUCE OR DISCLOSE ITS CONTENTS, SELL, LEASE, OR OTHERWISE
@rem TRANSFER IT TO ANY THIRD PARTY, IN WHOLE OR IN PART, WITHOUT THE SPECIFIC
@rem WRITTEN CONSENT OF THE MSC.SOFTWARE CORPORATION.
@rem
@setlocal
@echo off
set outnt=./plate_fem
set outnt=%outnt:/=\%
set archdir=c:/msc.software/msc_nastran_student_edition/2014/nastran/msc20140/win64i4
set archdir=%archdir:/=\%
set datecmd=c:/msc.software/msc_nastran_student_edition/2014/nastran/msc20140/win64i4/mscdate
set datecmd=%datecmd:/=\%
set exedir=
set jobnt=./plate_fem.T13504_52
set jobnt=%jobnt:/=\%
set lognt=./plate_fem.log
set lognt=%lognt:/=\%
set modedir=
set nastcmd=c:/msc.software/msc_nastran_student_edition/2014/nastran/msc20140/win64i4/nastran.exe
set nastcmd=%nastcmd:/=\%
set sdirnt=c:/scratch
set sdirnt=%sdirnt:/=\%
set solver=c:/msc.software/msc_nastran_student_edition/2014/nastran/msc20140/win64i4/analysis.exe
set solver=%solver:/=\%
set tcmdopt=
set tcmd=c:/msc.software/msc_nastran_student_edition/2014/nastran/msc20140/win64i4/msctime
set tcmd=%tcmd:/=\%
set user=Unknown
set MSC_BASE=c:/msc.software/msc_nastran_student_edition/2014/nastran
set MSC_VERSD=msc20140
set MSC_ARCH=win64i4
set MSC_BASE=%MSC_BASE:/=\%
set MSC_VERSD=%MSC_VERSD:/=\%
set MSC_ARCH=%MSC_ARCH:/=\%
set MSC_LICENSE_FILE=EDU
set MSC_LICENSE_FILE=%MSC_LICENSE_FILE:/=\%
set DBSDIR=c:/scratch
set JIDDIR=c:/users/hensomc/documents/matlab/vs_lam/vs_lam_working
set MSC_APP=no
set MSC_AUTH=EDU
set MSC_BASE=c:/msc.software/msc_nastran_student_edition/2014/nastran
set MSC_DBS=c:/scratch/plate_fem.T13504_52
set MSC_EXE=c:/msc.software/msc_nastran_student_edition/2014/nastran/msc20140/win64i4/analysis.exe
set MSC_ISHELLPATH=c:/users/hensomc/documents/matlab/vs_lam/vs_lam_working
set MSC_ISHELLEXT=
set MSC_JID=C:\Users\hensomc\Documents\MATLAB\vs_lam\vs_lam_working\plate_fem.bdf
set MSC_JIDPATH=
set MSC_JIDTYPE=.bdf
set MSC_MEM=26214400
set MSC_OLD=no
set MSC_OUT=./plate_fem
set MSC_SCR=yes
set MSC_SDIR=c:/scratch
set MSC_VERSD=msc20140
set OUTDIR=.
set TMP=c:/users/hensomc/appdata/local/temp
set MSC_SQAS_MESSAGE_FILE=c:/msc.software/msc_nastran_student_edition/2014/nastran\msc20140\win64i4\analysis.msg
if not "%DBSDIR%" == "" set DBSDIR=%DBSDIR:/=\%
if not "%JIDDIR%" == "" set JIDDIR=%JIDDIR:/=\%
if not "%MSC_BASE%" == "" set MSC_BASE=%MSC_BASE:/=\%
if not "%MSC_DBS%" == "" set MSC_DBS=%MSC_DBS:/=\%
if not "%MSC_EXE%" == "" set MSC_EXE=%MSC_EXE:/=\%
if not "%MSC_ISHELLPATH%" == "" set MSC_ISHELLPATH=%MSC_ISHELLPATH:/=\%
if not "%MSC_JID%" == "" set MSC_JID=%MSC_JID:/=\%
if not "%MSC_JIDPATH%" == "" set MSC_JIDPATH=%MSC_JIDPATH:/=\%
if not "%MSC_OUT%" == "" set MSC_OUT=%MSC_OUT:/=\%
if not "%MSC_SDIR%" == "" set MSC_SDIR=%MSC_SDIR:/=\%
if not "%OUTDIR%" == "" set OUTDIR=%OUTDIR:/=\%
if not "%TMP%" == "" set TMP=%TMP:/=\%
set MSC_ARCHD=%MSC_ARCH%
set MSC_ARCHM=%MSC_ARCH%
set HDF5_DISABLE_VERSION_CHECK=1
set MKL_PARDISO_OOC_CFG_PATH=%MSC_SDIR%
set MKL_PARDISO_OOC_CFG_FILE_NAME=
for %%f in (%MSC_DBS%) do set MKL_PARDISO_OOC_CFG_FILE_NAME=%%~nxf
set MKL_PARDISO_OOC_PATH=%MSC_SDIR%\%MKL_PARDISO_OOC_CFG_FILE_NAME%_ooc
set MKL_PARDISO_OOC_CFG_FILE_NAME=%MKL_PARDISO_OOC_CFG_FILE_NAME%_ooc.cfg
if "%MSC_SQAS_MESSAGE_FILE%" == "" (
   set MSC_SQAS_MESSAGE_FILE=%MSC_BASE%\%MSC_VERSD%\%MSC_ARCH%\analysis.msg
)
if "%SCA_SERVICE_CATALOG%" == "" (
  set SCA_SERVICE_CATALOG=%MSC_BASE%\%MSC_VERSD%\res\SCAServiceCatalog.xml
) else (
  set SCA_SERVICE_CATALOG=%SCA_SERVICE_CATALOG%;%MSC_BASE%\%MSC_VERSD%\res\SCAServiceCatalog.xml
)
if "%SCA_RESOURCE_DIR%" == "" (
  set SCA_RESOURCE_DIR=%MSC_BASE%\%MSC_VERSD%\res
) else (
  set SCA_RESOURCE_DIR=%SCA_RESOURCE_DIR%;%MSC_BASE%\%MSC_VERSD%\res
)
set dashes=-----------------------------------------------------------------
set findstr="%windir%\system32\findstr.exe"
set stars= ****************
set MSC_MSG=stderr
set MSC_SRV=
for /f "tokens=*" %%g in ('cd') do set CPWD=%%g
set PWD=C:\Users\hensomc\Documents\MATLAB\vs_lam\vs_lam_working
set PWD=%PWD:/=\%
cd /d %PWD%
for /f "tokens=*" %%g in ('"%datecmd%"') do set startdate=%%g
@rem
@rem Delete the old output files.
@rem
for %%f in (f06 f04 log ndb op2 pch plt asm becho marc.dat marc.log marc.out marc.sts marc.t16 marc.t19 des mnf PCS aeso sts out) do (
  if exist "%outnt%.%%f" del /f "%outnt%.%%f"
)
for %%f in (op2) do (
  setlocal EnableDelayedExpansion
  for %%i in (%outnt%) do set fname=%%~ni
  for %%i in (%outnt%) do set dpname=%%~dpi
  for /f "tokens=*" %%g in ('^( dir /b "%outnt%-*.%%f" 2^>nul: ^| %findstr% /r /i "!fname!-[0-9][0-9]*\.%%f" ^)') do (
    if exist "!dpname!%%g" del /f "!dpname!%%g" 2>nul:
  )  
  endlocal
)
for %%f in (log) do (
  if exist "%outnt%.digimat.%%f" del /f "%outnt%.digimat.%%f"
)
set asgf=%jobnt%.rcf
set PRINTLOC=^1^>^>"%lognt%" ^2^>^&^1
set PATH=%archdir%;%PATH%
if exist "%archdir%\lib\*.dll" set PATH=%archdir%\lib;%PATH%
if exist "%archdir%\nvidia\*.dll" set PATH=%archdir%\nvidia;%PATH%
if exist "%archdir%\nCode\bin\*.dll" set PYTHONPATH=%PYTHONPATH%;%archdir%\nCode\bin
if exist "%archdir%\nCode\bin\*.dll" set PATH=%archdir%\nCode\bin;%PATH%
if exist "%archdir%\lib\win32\*.dll" set PATH=%archdir%\lib\win32;%PATH%
if exist "%archdir%\lib\win64i4\*.dll" set PATH=%archdir%\lib\win64i4;%PATH%
if not "%SCA_LIBRARY_PATH%" == "" set PATH=%SCA_LIBRARY_PATH%;%PATH%
echo MSC Nastran V2014.0 (Intel Windows 7 Professional 6.1 7601) Control File: %PRINTLOC%
echo %dashes% %PRINTLOC%
type "%asgf%" %PRINTLOC%
echo %dashes% %PRINTLOC%
echo MSC_ISHELLPATH=%MSC_ISHELLPATH% %PRINTLOC%
echo MSC_JIDPATH=%MSC_JIDPATH% %PRINTLOC%
echo PATH=%PATH% %PRINTLOC%
if not "%SCA_LIBRARY_PATH%" == "" (
  echo SCA_LIBRARY_PATH=%SCA_LIBRARY_PATH% %PRINTLOC%
)
echo SCA_SERVICE_CATALOG=%SCA_SERVICE_CATALOG% %PRINTLOC%
echo SCA_RESOURCE_DIR=%SCA_RESOURCE_DIR% %PRINTLOC%
set SOL700=no
set mksol700=no
if "%SOL700%" == "yes" (
   if not exist "plate_fem.sol700.pth" (
      set mksol700=yes
      if "%MSC_ARCH%" == "win32" (
         echo %MSC_BASE%/%MSC_VERSD%/dyna/win32/run_dytran > plate_fem.sol700.pth
      ) else (
         echo %MSC_BASE%/%MSC_VERSD%/dyna/win64/run_dytran > plate_fem.sol700.pth
      )
      set  LSTC_VER=""
   )
)
set cmdline=C:\MSC.Software\MSC_Nastran_Student_Edition\2014\Nastran\msc20140\win64i4\nastran.exe plate_fem.bdf old=no scr=yes
echo Cmd Line:    %cmdline% %PRINTLOC%
echo Current Dir: %PWD% %PRINTLOC%
echo Time Cmd:    %tcmd%%tcmdopt% %PRINTLOC%
echo Executable:  %solver% %PRINTLOC%
if not "%MSC_USER_ARCH%" == "" echo Arch= Kywd:  %MSC_USER_ARCH%  %PRINTLOC%
if "%MSC_USER_ARCH%" == "" echo Arch= Kywd:  [Not specified]  %PRINTLOC%
echo MSC_MSG:     stderr  %PRINTLOC%
echo FPE:         yes  %PRINTLOC%
echo %dashes% %PRINTLOC%
set modeinfo=
if exist "%nastcmd%" (
   "%nastcmd%" -a0 -c0 limits %modeinfo% %PRINTLOC%
   echo %dashes% %PRINTLOC%
   "%nastcmd%" -a0 -c0 system %modeinfo% %PRINTLOC%
) else (
   echo Resource information data not available. %PRINTLOC%
)
echo %dashes% %PRINTLOC%
if not exist "%solver%" (
   echo Executable "%solver%" not found or not executable. %PRINTLOC%
   goto :jobdone
)
echo MSC Nastran started  %startdate% %PRINTLOC%
echo MSC Nastran started  %startdate%
set MSC_ASG=%asgf%
set MSC_MEM=26214400
"%tcmd%"%tcmdopt%    "%solver%"  FPE=yes   %PRINTLOC%
@rem
set /a f98 = 0
call :get_file_size  "%sdirnt%\plate_fem.T13504_52.f98"
if /i %FILE_SIZE% gtr 0 set /a f98 = 1
set /a f99 = 0
call :get_file_size  "%sdirnt%\plate_fem.T13504_52.f99"
if /i %FILE_SIZE% gtr 0 set /a f99 = 1
set /a i = %f98% + %f99%
if /i %i% neq 0 (
   for /F "tokens=*" %%f in ('^(type "%lognt%" ^| %findstr% /L NSEXIT ^)') do set fexit=%%f
   if not "%fexit%" == "" for /f "tokens=*" %%f in ('^( echo "%fexit%" ^| %findstr% /L "(0)" ^)') do set %ferr%=%%f
   if "%ferr%" == "" (
      if /i %f98% neq 0 (
         if /i %f99% neq 0 copy /f /y "%sdirnt%\plate_fem.T13504_52.f99"+"%sdirnt%\plate_fem.T13504_52.f98" "%outnt%.acms_out"
         if /i %f99% equ 0 copy /f /y "%sdirnt%\plate_fem.T13504_52.f98" "%outnt%.acms_out"
      ) else (
         copy /f /y "%sdirnt%\plate_fem.T13504_52.f99" "%outnt%.acms_out"
      )
   )
)
@rem
for %%f in ("%sdirnt%\plate_fem*_ts_*.dac") do call :delete_job_file "%%f"
for %%f in (dtout msg config) do call :concat_job_file %%f
for %%f in ("%sdirnt%\plate_fem.T13504_52.*") do call :delete_job_file "%%f"
for /f "tokens=*" %%g in ('"%datecmd%"') do set enddate=%%g
echo MSC Nastran finished %enddate% %PRINTLOC%
echo MSC Nastran finished %enddate%
for %%f in (ndb op2 pch plt xdb asm becho marc.dat marc.log marc.out marc.sts marc.t16 marc.t19 des mnf PCS aeso sts) do call :remove_0_file %%f
if exist "%MKL_PARDISO_OOC_CFG_FILE_NAME%" del /f "%MKL_PARDISO_OOC_CFG_FILE_NAME%"
if "%mksol700%" == "yes" (
   if exist plate_fem.sol700.pth del /f plate_fem.sol700.pth
)
@rem Delete job utility files
for %%f in ("%jobnt%*") do call :delete_job_file "%%f"
:jobdone
endlocal
goto :eof
@rem
@rem //////////////////////////////////////////////////////////////////
@rem
@rem   Function:  concat_job_file
@rem
:concat_job_file
setlocal
set ext=%~x1
type "plate_fem_*.%1" > plate_fem.%1 2>nul
del /f "plate_fem_*.%1" 2>nul:
>nul findstr "^" "plate_fem.%1" || del "plate_fem.%1"
endlocal
goto :eof
@rem
@rem
@rem //////////////////////////////////////////////////////////////////
@rem
@rem   Function:  delete_job_file
@rem
:delete_job_file
setlocal
set ext=%~x1
if "%ext%" == ".bat" goto :delete_job_file_done
if "%ext%" == ".rcf" goto :delete_job_file_done
del /f "%1" 2>nul:
:delete_job_file_done
endlocal
goto :eof
@rem
@rem
@rem //////////////////////////////////////////////////////////////////
@rem
@rem   Function:  get_file_size
@rem
:get_file_size
setlocal
set /a i = -1
if not %1 == "" for /f "tokens=4" %%f in ('^dir /-c %1 2^>nul:^)') do set i="%%f"
endlocal && set FILE_SIZE=%i%
goto :eof
@rem
@rem //////////////////////////////////////////////////////////////////
@rem
@rem   Function:  remove_0_file
@rem
:remove_0_file
if not exist "%outnt%.%1" goto :remove_0_file_done
call :get_file_size "%outnt%.%1"
if /i "%FILE_SIZE%" equ 0 del /f "%outnt%.%1" 2>nul:
:remove_0_file_done
goto :eof
@rem
@rem //////////////////////////////////////////////////////////////////
