rem update file src\General\svn_rev.f90 with current revision numbering to be included in the binary
rem to be executed within the directory where you downloaded WASA-SED files into (manually or by Makefile)

@echo off
rem echo "called with %1"
set src_path=%1
if "%src_path%" == "" set src_path=.

echo !this file is updated by update_revision_no.sh/.bat > %src_path%\General\svn_rev.f90

rem set output of git to variable
for /f %%i in ('git describe --tags --long') do set gitstr=%%i 

rem check for local modifications
set cmmt=
for /f %%a in ('git status --porcelain') do set mod_files="%%a"
if "%mod_files%" NEQ "" set cmmt=(local modifications)

echo rev_string1='%gitstr% %cmmt%' >>  %src_path%\General\svn_rev.f90

rem set output of git to variable
rem for /f %%i in ('git show -s --format^=%%ci HEAD') do set gitdate=%%i 
rem for /f "delims=" %%a in ('git show -s --format^=%%ci %1^^') do set hours=!hours! "%%a"
for /f "delims=" %%a in ('git show -s --format^=%%ci ^%1^') do set gitdate=%%a
rem for /f "delims=" %%a in ('git show -s --format^=%%ci ^%1') do set gitdate=%%a

rem for /f "delims=" %a in ('git show -s --format ^=%ci') do set hours=!hours! "%a"

for /f %%i in ('date /T') do set sysdate=%%i 
for /f %%i in ('time /T') do set systime=%%i 

echo rev_string2='repository date %gitdate%, built %sysdate%%systime%' >>  %src_path%\General\svn_rev.f90


