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
# ADD F90 /browser /check:bounds /check:power /check:overflow /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:declarations /warn:nofileopt
# SUBTRACT F90 /check:underflow /fltconsistency /threads
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
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

SOURCE=.\General\general_h.f90
# End Source File
# End Group
# Begin Source File

SOURCE=.\General\calcday.f90
DEP_F90_CALCD=\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\General\calcyear.f90
DEP_F90_CALCY=\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\General\climo.f90
DEP_F90_CLIMO=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\time_h.mod"\
	".\Debug\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\General\petcalc.f90
DEP_F90_PETCA=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\General\readgen.f90
DEP_F90_READG=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	".\General\allocat_general.var"\
	".\Hillslope\allocat_erosion.var"\
	".\Hillslope\allocat_hymo.var"\
	".\Reservoir\allocat_reservoir_lake.var"\
	".\River\allocat_routing.var"\
	
# End Source File
# Begin Source File

SOURCE=.\General\utils_h.f90
# End Source File
# Begin Source File

SOURCE=.\General\wasa.f90
DEP_F90_WASA_=\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
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
# End Source File
# End Group
# Begin Source File

SOURCE=.\Hillslope\check_climate.f90
DEP_F90_CHECK=\
	".\Debug\climo_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\etp_max.f90
DEP_F90_ETP_M=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\etp_soil.f90
DEP_F90_ETP_S=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\etp_soil_hour.f90
DEP_F90_ETP_SO=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\hymo_all.f90
DEP_F90_HYMO_=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\model_state_io.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	".\Debug\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\model_state_io.f90
DEP_F90_MODEL=\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\readhymo.f90
DEP_F90_READH=\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	".\Debug\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\sedi_yield.f90
DEP_F90_SEDI_=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\sedi_yield_subbas.f90
DEP_F90_SEDI_Y=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\soildistr.f90
DEP_F90_SOILD=\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hillslope\soilwat.f90
DEP_F90_SOILW=\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\utils_h.mod"\
	
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
	".\Debug\common_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\muskingum.f90
DEP_F90_MUSKI=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\route_sediments.f90
DEP_F90_ROUTE=\
	".\Debug\common_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\routing.f90
DEP_F90_ROUTI=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\routing_coefficients.f90
DEP_F90_ROUTIN=\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\routing_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\River\routing_new.f90
DEP_F90_ROUTING=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
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

SOURCE=.\Reservoir\reservoir_lake_h.f90
# End Source File
# End Group
# Begin Source File

SOURCE=.\Reservoir\change_sec.f90
DEP_F90_CHANG=\
	".\Debug\common_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq1_wu.f90
DEP_F90_EQ1_W=\
	".\Debug\common_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq2_ashida.f90
DEP_F90_EQ2_A=\
	".\Debug\common_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq3_tsinghua.f90
DEP_F90_EQ3_T=\
	".\Debug\common_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\eq4_ackers.f90
DEP_F90_EQ4_A=\
	".\Debug\common_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\hydraul_res.f90
DEP_F90_HYDRA=\
	".\Debug\common_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\lake.f90
DEP_F90_LAKE_=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\erosion_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	".\Debug\utils_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\lake_routing.f90
DEP_F90_LAKE_R=\
	".\Debug\common_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\reservoir.f90
DEP_F90_RESER=\
	".\Debug\climo_h.mod"\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\reservoir_routing.f90
DEP_F90_RESERV=\
	".\Debug\common_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\sedbal.f90
DEP_F90_SEDBA=\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\sedbal_lake.f90
DEP_F90_SEDBAL=\
	".\Debug\common_h.mod"\
	".\Debug\lake_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\semres.f90
DEP_F90_SEMRE=\
	".\Debug\common_h.mod"\
	".\Debug\hymo_h.mod"\
	".\Debug\params_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\routing_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Reservoir\vert_dist.f90
DEP_F90_VERT_=\
	".\Debug\common_h.mod"\
	".\Debug\reservoir_h.mod"\
	".\Debug\time_h.mod"\
	
# End Source File
# End Group
# Begin Group "Parameter"

# PROP Default_Filter ""
# End Group
# End Target
# End Project
