# Microsoft Developer Studio Project File - Name="WASA" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=WASA - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "WASA.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "WASA.mak" CFG="WASA - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "WASA - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "WASA - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "WASA - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# SUBTRACT F90 /fast
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "WASA - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir ""
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /check:format /check:output_conversion /compile_only /dbglibs /debug:full /fpscomp:general /nologo /reentrancy:threaded /stand:none /traceback /warn:argument_checking /warn:declarations /warn:nofileopt /warn:truncated_source /warn:unused
# SUBTRACT F90 /check:arg_temp_created
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /profile /debug /machine:I386

!ENDIF 

# Begin Target

# Name "WASA - Win32 Release"
# Name "WASA - Win32 Debug"
# Begin Group "General"

# PROP Default_Filter ""
# Begin Group "General_Header"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\General\climo_h.f90
# End Source File
# Begin Source File

SOURCE=.\General\common_h.f90
# End Source File
# Begin Source File

SOURCE=.\General\params_h.f90
# End Source File
# Begin Source File

SOURCE=.\General\svn_rev.var
# End Source File
# Begin Source File

SOURCE=.\General\time_h.f90
# End Source File
# Begin Source File

SOURCE=.\General\utils_h.f90
# End Source File
# End Group
# Begin Source File

SOURCE=.\General\calcyear.f90
DEP_F90_CALCY=\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\General\climo.f90
DEP_F90_CLIMO=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\General\petcalc.f90
DEP_F90_PETCA=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\params_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\General\readgen.f90
DEP_F90_READG=\
	".\General\allocat_general.var"\
	".\Hillslope\allocat_erosion.var"\
	".\Hillslope\allocat_hymo.var"\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	".\Reservoir\allocat_reservoir_lake.var"\
	".\River\allocat_routing.var"\
	
# End Source File
# Begin Source File

SOURCE=.\General\wasa.f90
DEP_F90_WASA_=\
	".\General\svn_rev.var"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# End Group
# Begin Group "Hillslope"

# PROP Default_Filter ""
# Begin Group "Hillslope_Header"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Hillslope\allocat_erosion.var
# End Source File
# Begin Source File

SOURCE=.\Hillslope\allocat_hymo.var
# End Source File
# Begin Source File

SOURCE=.\Hillslope\erosion_h.f90
# End Source File
# Begin Source File

SOURCE=.\Hillslope\hymo_h.f90
DEP_F90_HYMO_=\
	".\Release\common_h.mod"\
	".\Release\params_h.mod"\
	
# End Source File
# End Group
# Begin Source File

SOURCE=.\Hillslope\check_climate.f90
DEP_F90_CHECK=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\etp_max.f90
DEP_F90_ETP_M=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\etp_soil.f90
DEP_F90_ETP_S=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\etp_soil_hour.f90
DEP_F90_ETP_SO=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\hymo_all.f90
DEP_F90_HYMO_A=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\model_state_io.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\model_state_io.f90
DEP_F90_MODEL=\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\readhymo.f90
DEP_F90_READH=\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\sedi_yield.f90
DEP_F90_SEDI_=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\sedi_yield_subbas.f90
DEP_F90_SEDI_Y=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\soildistr.f90
DEP_F90_SOILD=\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\soilwat.f90
DEP_F90_SOILW=\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# End Group
# Begin Group "River"

# PROP Default_Filter ""
# Begin Group "River_Header"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\River\allocat_routing.var
# End Source File
# Begin Source File

SOURCE=.\River\routing_h.f90
# End Source File
# End Group
# Begin Source File

SOURCE=.\River\bedload.f90
DEP_F90_BEDLO=\
	".\Release\common_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\muskingum.f90
DEP_F90_MUSKI=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\route_sediments.f90
DEP_F90_ROUTE=\
	".\Release\common_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\routing.f90
DEP_F90_ROUTI=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\routing_coefficients.f90
DEP_F90_ROUTIN=\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\model_state_io.mod"\
	".\Release\routing_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\routing_new.f90
DEP_F90_ROUTING=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# End Group
# Begin Group "Reservoir"

# PROP Default_Filter ""
# Begin Group "Reservoir_Header"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Reservoir\allocat_reservoir_lake.var
# End Source File
# Begin Source File

SOURCE=.\Reservoir\lake_h.f90
# End Source File
# Begin Source File

SOURCE=.\Reservoir\reservoir_h.f90
# End Source File
# End Group
# Begin Source File

SOURCE=.\Reservoir\change_sec.f90
DEP_F90_CHANG=\
	".\Release\common_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq1_wu.f90
DEP_F90_EQ1_W=\
	".\Release\common_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq2_ashida.f90
DEP_F90_EQ2_A=\
	".\Release\common_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq3_tsinghua.f90
DEP_F90_EQ3_T=\
	".\Release\common_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq4_ackers.f90
DEP_F90_EQ4_A=\
	".\Release\common_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\hydraul_res.f90
DEP_F90_HYDRA=\
	".\Release\common_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\lake.f90
DEP_F90_LAKE_=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\erosion_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\lake_routing.f90
DEP_F90_LAKE_R=\
	".\Release\common_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\reservoir.f90
DEP_F90_RESER=\
	".\Release\climo_h.mod"\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	".\Release\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\reservoir_routing.f90
DEP_F90_RESERV=\
	".\Release\common_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\sedbal.f90
DEP_F90_SEDBA=\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\sedbal_lake.f90
DEP_F90_SEDBAL=\
	".\Release\common_h.mod"\
	".\Release\lake_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\semres.f90
DEP_F90_SEMRE=\
	".\Release\common_h.mod"\
	".\Release\hymo_h.mod"\
	".\Release\params_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\routing_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\vert_dist.f90
DEP_F90_VERT_=\
	".\Release\common_h.mod"\
	".\Release\reservoir_h.mod"\
	".\Release\time_h.mod"\
	
# End Source File
# End Group
# Begin Group "Parameter"

# PROP Default_Filter ""
# End Group
# End Target
# End Project
