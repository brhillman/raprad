F77 = f77

FFLAGS = -g -c

LINK.libs= -lc -lF77 -lV77 -lm -lM77 -lsunmath



OBJECTS1=\
	mk_sunview_st_atm.o




mk_sunview_st_atm: ${OBJECTS1}
	$(F77) -g -o $@ ${OBJECTS1} ${LINK.libs}





mk_sunview_st_atm.o:
	$(F77) ${FFLAGS} mk_sunview_st_atm.f


