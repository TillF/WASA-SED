Tobias Pilz <tpilz@uni-potsdam.de>, 28. July 2014

compiled code with gfortran 4.8 with flags: (see Makefile or compiler documentation for descriptions)
-Wno-maybe-uninitialized
-Wtabs
-fbacktrace
-g
-fcheck=all
-Wall
-Wextra
-fimplicit-none
-ffree-line-length-none

warnings not yet fixed:
- lots of warnings like: 'Warning: Equality comparison for REAL(4) at (1)'
	- not sure if this is relevant
	- quite new to gfortran (since 2012?) within -Wall flag
	- in internet: discussions whether this is important but it might cause unexpected behaviour
	- unexperienced programmers should resolve the warnings
- some warnings like: 'Warning: Legacy Extension: Label at (1) is not in the same block as the GOTO statement at (2)'
	- not sure if this is relevant
- hymo_all.f90 line 975: runtime warning: An array temporary was created
	- not sure what to do and what that means
