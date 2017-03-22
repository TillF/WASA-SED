rem SubWCRev .\ General\svn_rev_template.txt General\svn_rev.var
echo !this file is updated on calling update_revision_no.bat/.sh to create svn_rev.var which is linked to the source code > General\svn_rev.var

rem set output of git to variable
for /f %%i in ('git describe --tags --long') do set gitstr=%%i 

rem check for local modifications
cmmt=
for /f %%a in ('git status --porcelain') do set mod_files="%%a"
if "%mod_files%" NEQ "" set cmmt=(local modifications)

echo rev_string1='%gitstr% %cmmt%' >> General\svn_rev.var

rem set output of git to variable
rem for /f %%i in ('git show -s --format^=%%ci HEAD') do set gitdate=%%i 
rem for /f "delims=" %%a in ('git show -s --format^=%%ci %1^^') do set hours=!hours! "%%a"
for /f "delims=" %%a in ('git show -s --format^=%%ci ^%1^') do set gitdate=%%a

rem for /f "delims=" %a in ('git show -s --format ^=%ci') do set hours=!hours! "%a"

for /f %%i in ('date /T') do set sysdate=%%i 
for /f %%i in ('time /T') do set systime=%%i 

echo rev_string2='repository date %gitdate%, built %sysdate%%systime%' >> General\svn_rev.var


