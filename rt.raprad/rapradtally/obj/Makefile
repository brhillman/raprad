
#
#
#          Makefile generator for Raprad radiative transfer code.
#                   Seiji Kato 11 December 1996
#                     Penn State University
#
#


# SUN System compilier options
SUNCC = cc -g -I/usr/local/include -c -o $@

# GNU System compilier options
GNUCC = gcc -traceback -I/usr/local/include

# SGI System compilier options
SGICC = cc -g -cckr

#F77 = f77
F77 = ifort -traceback

FFLAGS = -g -c -o

#LINK.c= -lc -lF77 -lM77 -lV77 -lm -lsunmath
LINK.c= -lc -lm -nofor-main


COBJECTS = \
	 ../../util/obj/atmosphere.layers.interpolate.o \
	 ../../util/obj/atmosphere.layers.subdivide.o \
	 ../../util/obj/atmosphere.layers.subdivide.check.o \
	 ../../util/obj/atmosphere.read.mcclatchey.o \
	 ../../util/obj/cntnmmodel.check.o \
	 ../../util/obj/cntnmmodel.read.o \
	 ../../util/obj/numrec.double.piksrt.o \
	 ../../util/obj/numrec.nrutil.o \
	 ../../util/obj/parsers.parameters.o \
	 ../../util/obj/raprad.check.constituents.o \
	 ../../util/obj/raprad.check.geometry.sunzenith.o \
	 ../../util/obj/raprad.check.photon.partition.o \
	 ../../util/obj/raprad.check.photon.space.final.o \
	 ../../util/obj/raprad.photon.space.final.o \
	 ../../util/obj/raprad.read.configuration.o \
	 ../../util/obj/raprad.read.constituents.o \
	 ../../util/obj/raprad.read.geometry.sunzenith.o \
	 ../../util/obj/raprad.read.photon.partition.o \
	 ../../util/obj/raprad.read.photon.space.o \
	 ../../util/obj/spectralmodel.check.kato.o \
	 ../../util/obj/spectralmodel.read.kato.o \
	 ../../util/obj/spectralmodel.check.pollack.o \
	 ../../util/obj/spectralmodel.read.pollack.o \
	 ../../util/obj/spectralmodel.check.mlawer_lw.o \
	 ../../util/obj/spectralmodel.read.mlawer_lw.o \
	 ../../util/obj/thermodynamics.units.o \
	 raprad.bbflux.main.o



all:			\
	rapradbbfluxlw 	\
	rapradbbfluxsw  \
        combine_sw_lw_flx
        


OBJECTS2 = \
	combine_sw_lw_flx.o



rapradbbfluxlw: raprad.bbflux.lw.o $(COBJECTS)
	 $(F77) -g -o $@ raprad.bbflux.lw.o $(COBJECTS) $(LINK.c)

rapradbbfluxsw: raprad.bbflux.sw.o $(COBJECTS)
	 $(F77) -g -o $@ raprad.bbflux.sw.o $(COBJECTS) $(LINK.c)

read_merged_sounding: ${OBJECTS1}
	$(F77) -g -o $@ ${OBJECTS1} ${CDFLIB01} ${LINK.libs}

combine_sw_lw_flx: $(OBJECTS2)
	 $(F77) -g -o $@ combine_sw_lw_flx $(OBJECTS2) $(LINK.c)


raprad.bbflux.main.o: ../src/raprad.bbflux.main.c
	 $(SUNCC) ../src/raprad.bbflux.main.c

raprad.bbflux.lw.o: ../src/raprad.bbflux.lw.f
	$(F77) $(FFLAGS) $@ ../src/raprad.bbflux.lw.f

raprad.bbflux.sw.o: ../src/raprad.bbflux.sw.f
	$(F77) $(FFLAGS) $@ ../src/raprad.bbflux.sw.f


combine_sw_lw_flx.o: ../src/combine_sw_lw_flx.f
	$(F77) $(FFLAGS) $@ ../src/combine_sw_lw_flx.f


../../util/obj/raprad.check.constituents.o: ../../util/src/raprad.check.constituents.c
	 $(SUNCC) ../../util/src/raprad.check.constituents.c


../../util/obj/raprad.check.geometry.sunzenith.o: ../../util/src/raprad.check.geometry.sunzenith.c
	 $(SUNCC) ../../util/src/raprad.check.geometry.sunzenith.c


../../util/obj/raprad.check.photon.partition.o: ../../util/src/raprad.check.photon.partition.c
	 $(SUNCC) ../../util/src/raprad.check.photon.partition.c


../../util/obj/raprad.check.photon.space.final.o: ../../util/src/raprad.check.photon.space.final.c
	 $(SUNCC) ../../util/src/raprad.check.photon.space.final.c


../../util/obj/atmosphere.layers.interpolate.o: ../../util/src/atmosphere.layers.interpolate.c
	 $(SUNCC) ../../util/src/atmosphere.layers.interpolate.c


../../util/obj/atmosphere.layers.subdivide.o: ../../util/src/atmosphere.layers.subdivide.c
	 $(SUNCC) ../../util/src/atmosphere.layers.subdivide.c


../../util/obj/atmosphere.layers.subdivide.check.o: ../../util/src/atmosphere.layers.subdivide.check.c
	 $(SUNCC) ../../util/src/atmosphere.layers.subdivide.check.c


../../util/obj/numrec.double.piksrt.o: ../../util/src/numrec.double.piksrt.c
	 $(SUNCC) ../../util/src/numrec.double.piksrt.c


../../util/obj/numrec.nrutil.o: ../../util/src/numrec.nrutil.c
	 $(SUNCC) ../../util/src/numrec.nrutil.c


../../util/obj/parsers.parameters.o: ../../util/src/parsers.parameters.c
	 $(SUNCC) ../../util/src/parsers.parameters.c


../../util/obj/raprad.photon.space.final.o: ../../util/src/raprad.photon.space.final.c
	 $(SUNCC) ../../util/src/raprad.photon.space.final.c


../../util/obj/raprad.read.configuration.o: ../../util/src/raprad.read.configuration.c
	 $(SUNCC) ../../util/src/raprad.read.configuration.c


../../util/obj/raprad.read.constituents.o: ../../util/src/raprad.read.constituents.c
	 $(SUNCC) ../../util/src/raprad.read.constituents.c


../../util/obj/raprad.read.geometry.sunzenith.o: ../../util/src/raprad.read.geometry.sunzenith.c
	 $(SUNCC) ../../util/src/raprad.read.geometry.sunzenith.c


../../util/obj/atmosphere.read.mcclatchey.o: ../../util/src/atmosphere.read.mcclatchey.c
	 $(SUNCC) ../../util/src/atmosphere.read.mcclatchey.c


../../util/obj/cntnmmodel.check.o: ../../util/src/cntnmmodel.check.c
	 $(SUNCC) ../../util/src/cntnmmodel.check.c


../../util/obj/cntnmmodel.read.o: ../../util/src/cntnmmodel.read.c
	 $(SUNCC) ../../util/src/cntnmmodel.read.c


../../util/obj/raprad.read.photon.partition.o: ../../util/src/raprad.read.photon.partition.c
	 $(SUNCC) ../../util/src/raprad.read.photon.partition.c


../../util/obj/raprad.read.photon.space.o: ../../util/src/raprad.read.photon.space.c
	 $(SUNCC) ../../util/src/raprad.read.photon.space.c


../../util/obj/spectralmodel.check.kato.o: ../../util/src/spectralmodel.check.kato.c
	 $(SUNCC) ../../util/src/spectralmodel.check.kato.c


../../util/obj/spectralmodel.check.pollack.o: ../../util/src/spectralmodel.check.pollack.c
	 $(SUNCC) ../../util/src/spectralmodel.check.pollack.c


../../util/obj/spectralmodel.read.kato.o: ../../util/src/spectralmodel.read.kato.c
	 $(SUNCC) ../../util/src/spectralmodel.read.kato.c


../../util/obj/spectralmodel.read.pollack.o: ../../util/src/spectralmodel.read.pollack.c
	 $(SUNCC) ../../util/src/spectralmodel.read.pollack.c


../../util/obj/spectralmodel.read.mlawer_lw.o: ../../util/src/spectralmodel.read.mlawer_lw.c
	 $(SUNCC) ../../util/src/spectralmodel.read.mlawer_lw.c


../../util/obj/spectralmodel.check.mlawer_lw.o: ../../util/src/spectralmodel.check.mlawer_lw.c
	 $(SUNCC) ../../util/src/spectralmodel.check.mlawer_lw.c


../../util/obj/thermodynamics.units.o: ../../util/src/thermodynamics.units.c
	 $(SUNCC) ../../util/src/thermodynamics.units.c



clean:
	 rm $(COBJECTS) $(OBJECTS1) $(OBJECTS2)
	 rm raprad.bbflux.lw.o
	 rm raprad.bbflux.sw.o



#
#          END OF MAKE
#

