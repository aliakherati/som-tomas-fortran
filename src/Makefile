# Makefile for TOMAS box model
# questions or comments, email Jeff Pierce, jrpierce@andrew.cmu.edu
###########################################################################

# options
comp = ifort
# to run quickly, use -O2 flag (optimize level 2) or -O if using g77 with no -g flag (for debugger)
# to debug, remove -O2 and add -g flag
o_flags = -c -O3 # flags for the compiling of *.o files
exe_flags = -O3 #-static # flags for the compiling of *.exe files


# groups
box_obj  = box.o	boxinit.o	getlen.o	wallloss.o
cond_obj = tmcond.o 	dmdt_int.o  ezcond.o    ezmass.o    \
           getCondSink.o      getMassSink.o     soacond.o
coag_obj = multicoag.o
size_obj = initbounds.o	loginit.o	report.o	\
	   storenm.o	ezwatereqm.o	eznh3eqm.o	\
	   mnfix.o	waterso4.o	waternacl.o
prop_obj = gasdiff.o	aerodens.o
som_obj  = saprc.o	integr2.o	rhs.o	saprc14_rev1.o

# executables
box.exe : $(box_obj) $(cond_obj) $(coag_obj) $(size_obj) $(prop_obj) $(som_obj) sizecode.COM modlspc.h gaskin.h common.inc Makefile
	$(comp)  $(exe_flags) -o box.exe $(box_obj) $(cond_obj) $(coag_obj) $(size_obj) $(prop_obj) $(som_obj)

# object files
box.o : box.f sizecode.COM modlspc.h gaskin.h common.inc Makefile
	$(comp) $(o_flags) box.f
wallloss.o: wallloss.f sizecode.COM Makefile
	$(comp) $(o_flags) wallloss.f
boxinit.o : boxinit.f sizecode.COM Makefile
	$(comp) $(o_flags) boxinit.f
tmcond.o : tmcond.f sizecode.COM Makefile
	$(comp) $(o_flags) tmcond.f
ezcond.o : ezcond.f sizecode.COM Makefile
	$(comp) $(o_flags) ezcond.f
ezmass.o : ezmass.f sizecode.COM Makefile
	$(comp) $(o_flags) ezmass.f
getCondSink.o : getCondSink.f sizecode.COM Makefile
	$(comp) $(o_flags) getCondSink.f
getMassSink.o : getMassSink.f sizecode.COM Makefile
	$(comp) $(o_flags) getMassSink.f
getlen.o : getlen.f Makefile
	$(comp) $(o_flags) getlen.f
multicoag.o : multicoag.f sizecode.COM Makefile
	$(comp) $(o_flags) multicoag.f
soacond.o : soacond.f sizecode.COM Makefile
	$(comp) $(o_flags) soacond.f
dmdt_int.o : dmdt_int.f Makefile
	$(comp) $(o_flags) dmdt_int.f
initbounds.o : initbounds.f sizecode.COM Makefile
	$(comp) $(o_flags) initbounds.f
loginit.o : loginit.f sizecode.COM Makefile
	$(comp) $(o_flags) loginit.f
report.o : report.f sizecode.COM Makefile
	$(comp) $(o_flags) report.f
storenm.o : storenm.f sizecode.COM Makefile
	$(comp) $(o_flags) storenm.f
ezwatereqm.o : ezwatereqm.f sizecode.COM Makefile
	$(comp) $(o_flags) ezwatereqm.f
eznh3eqm.o : eznh3eqm.f sizecode.COM Makefile
	$(comp) $(o_flags) eznh3eqm.f
gasdiff.o  : gasdiff.f Makefile
	$(comp) $(o_flags) gasdiff.f
aerodens.o : aerodens.f Makefile
	$(comp) $(o_flags) aerodens.f
mnfix.o : mnfix.f sizecode.COM Makefile
	$(comp) $(o_flags) mnfix.f
waterso4.o : waterso4.f Makefile
	$(comp) $(o_flags) waterso4.f
waternacl.o : waternacl.f Makefile
	$(comp) $(o_flags) waternacl.f
saprc.o : saprc.f modlspc.h Makefile
	$(comp) $(o_flags) saprc.f
integr2.o : integr2.f modlspc.h common.inc Makefile
	$(comp) $(o_flags) integr2.f
rhs.o : rhs.f modlspc.h gaskin.h common.inc Makefile
	$(comp) $(o_flags) rhs.f
saprc14_rev1.o : saprc14_rev1.f Makefile
	$(comp) $(o_flags) saprc14_rev1.f

clean :
	rm *.o box.exe
