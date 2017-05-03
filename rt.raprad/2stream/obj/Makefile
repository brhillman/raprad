
#
#
#          Makefile generator for Raprad radiative transfer code.
# modified by       Seiji Kato       11 December 1996
#                   Eugene Clothiaux 27 May      1998
#                     Penn State University
#
#


# SUN System compilier options
SUNCC = cc -g -I/usr/local/include -c -o $@

# GNU System compilier options
GNUCC = gcc -I/usr/local/include

# SGI System compilier options
SGICC = cc -g -cckr

#F77 = f77
F77 = ifort -traceback

FFLAGS = -g -c -o $@

#LINK.c= -lc -lF77 -lM77 -lV77 -lm -lsunmath
LINK.c= -lc -lm -nofor-main

OBJECTS = \
	 ../../util/obj/atmosphere.clausiusclapyron.o \
	 ../../util/obj/atmosphere.gasamounts.bext.o \
	 ../../util/obj/atmosphere.gasamounts.u.o \
	 ../../util/obj/atmosphere.gasamounts.u.check.o \
	 ../../util/obj/atmosphere.layers.interpolate.o \
	 ../../util/obj/atmosphere.layers.subdivide.o \
	 ../../util/obj/atmosphere.layers.subdivide.check.o \
	 ../../util/obj/atmosphere.mfp.aerosol.o \
	 ../../util/obj/atmosphere.mfp.cloud.o \
	 ../../util/obj/atmosphere.mfp.spectralmodel.kato.o \
	 ../../util/obj/atmosphere.mfp.spectralmodel.mlawer_lw.o \
	 ../../util/obj/rrtm_driver_setcoef_taumol.o \
	 ../../util/obj/setcoef.o \
	 ../../util/obj/taumol.o \
	 ../../util/obj/k_g.o \
	 ../../util/obj/atmosphere.mfp.spectralmodel.pollack.o \
	 ../../util/obj/atmosphere.mfp.ozone.martonchik.o \
	 ../../util/obj/atmosphere.mfp.rayleigh.o \
	 ../../util/obj/atmosphere.mfp.rayleigh.martonchik.o \
	 ../../util/obj/atmosphere.mfp.template.o \
	 ../../util/obj/atmosphere.read.mcclatchey.o \
	 ../../util/obj/cntnm.xsect.o \
	 ../../util/obj/cntnmmodel.check.o \
	 ../../util/obj/cntnmmodel.read.o \
	 ../../util/obj/numrec.double.piksrt.o \
	 ../../util/obj/numrec.nrutil.o \
	 ../../util/obj/parsers.parameters.o \
	 ../../util/obj/phasefcns.henyeygreenstein.o \
	 ../../util/obj/phasefcns.read.ascii.o \
	 ../../util/obj/phasefcns.read.binary.o \
	 ../../util/obj/raprad.alloc.o \
	 ../../util/obj/raprad.atmosphere.layers.o \
	 ../../util/obj/raprad.atmosphere.layers.tau.o \
	 ../../util/obj/raprad.check.atmosphere_layers.o \
	 ../../util/obj/raprad.check.atmosphere_layers.lw.o \
	 ../../util/obj/raprad.check.brdf.selection.o \
	 ../../util/obj/raprad.check.constituents.o \
	 ../../util/obj/raprad.check.geometry.sunzenith.o \
	 ../../util/obj/raprad.check.photon.partition.o \
	 ../../util/obj/raprad.check.photon.space.final.o \
	 ../../util/obj/raprad.photon.space.final.o \
	 ../../util/obj/raprad.read.brdf.selection.o \
	 ../../util/obj/raprad.read.configuration.o \
	 ../../util/obj/raprad.read.constituent.phasefcn.o \
	 ../../util/obj/raprad.read.constituents.o \
	 ../../util/obj/raprad.read.geometry.sunzenith.o \
	 ../../util/obj/raprad.read.photon.partition.o \
	 ../../util/obj/raprad.read.photon.space.o \
	 ../../util/obj/rayleigh.bext.o \
	 ../../util/obj/rayleigh.crosssection.o \
	 ../../util/obj/rayleigh.indexrefraction.air.o \
	 ../../util/obj/rayleigh.meanfreepath.o \
	 ../../util/obj/rayleigh.opticaldepth.o \
	 ../../util/obj/spectralmodel.check.kato.o \
	 ../../util/obj/spectralmodel.check.mlawer_lw.o \
	 ../../util/obj/spectralmodel.check.pollack.o \
	 ../../util/obj/spectralmodel.read.kato.o \
	 ../../util/obj/spectralmodel.read.mlawer_lw.o \
	 ../../util/obj/spectralmodel.read.pollack.o \
	 ../../util/obj/thermodynamics.clausiusclapyron.o \
	 ../../util/obj/thermodynamics.units.o


2STREAMLWOBJECTS = \
	raprad.2stream.lw.main.o \
	raprad.2stream.lw.add.o \
	raprad.2stream.lw.driver.o \
	raprad.2stream.lw.twostr.o \
	raprad.2stream.lw.oppr1.o \
	raprad.2stream.lw.newflux1.o 


2STREAMSWOBJECTS = \
	raprad.2stream.sw.main.o \
	raprad.2stream.sw.add.o \
	raprad.2stream.sw.driver.o \
	raprad.2stream.sw.twostr.o 


2STREAMLWMIXOBJECTS = \
	raprad.2stream.lw_mix.main.o \
	raprad.2stream.lw_mix.o  


all: \
	raprad2streamlw \
	raprad2streamsw \
#	raprad2streamlw_mix


raprad2streamlw: $(2STREAMLWOBJECTS) $(OBJECTS)
	$(F77) -g -o $@ $(2STREAMLWOBJECTS) $(OBJECTS) $(LINK.c)


raprad2streamsw: $(2STREAMSWOBJECTS) $(OBJECTS)
	$(F77) -g -o $@ $(2STREAMSWOBJECTS) $(OBJECTS) $(LINK.c)


raprad2streamlw_mix: $(2STREAMLWMIXOBJECTS) $(OBJECTS)
	$(F77) -g -o $@ $(2STREAMLWMIXOBJECTS) $(OBJECTS) $(LINK.c)


../../util/obj/atmosphere.clausiusclapyron.o: ../../util/src/atmosphere.clausiusclapyron.c
	 $(SUNCC) ../../util/src/atmosphere.clausiusclapyron.c


../../util/obj/atmosphere.gasamounts.bext.o: ../../util/src/atmosphere.gasamounts.bext.c
	 $(SUNCC) ../../util/src/atmosphere.gasamounts.bext.c


../../util/obj/atmosphere.gasamounts.u.o: ../../util/src/atmosphere.gasamounts.u.c
	 $(SUNCC) ../../util/src/atmosphere.gasamounts.u.c


../../util/obj/atmosphere.gasamounts.u.check.o: ../../util/src/atmosphere.gasamounts.u.check.c
	 $(SUNCC) ../../util/src/atmosphere.gasamounts.u.check.c


../../util/obj/atmosphere.layers.interpolate.o: ../../util/src/atmosphere.layers.interpolate.c
	 $(SUNCC) ../../util/src/atmosphere.layers.interpolate.c


../../util/obj/atmosphere.layers.subdivide.o: ../../util/src/atmosphere.layers.subdivide.c
	 $(SUNCC) ../../util/src/atmosphere.layers.subdivide.c


../../util/obj/atmosphere.layers.subdivide.check.o: ../../util/src/atmosphere.layers.subdivide.check.c
	 $(SUNCC) ../../util/src/atmosphere.layers.subdivide.check.c


../../util/obj/atmosphere.mfp.aerosol.o: ../../util/src/atmosphere.mfp.aerosol.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.aerosol.c


../../util/obj/atmosphere.mfp.cloud.o: ../../util/src/atmosphere.mfp.cloud.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.cloud.c


../../util/obj/atmosphere.mfp.ozone.martonchik.o: ../../util/src/atmosphere.mfp.ozone.martonchik.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.ozone.martonchik.c


../../util/obj/atmosphere.mfp.rayleigh.o: ../../util/src/atmosphere.mfp.rayleigh.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.rayleigh.c


../../util/obj/atmosphere.mfp.rayleigh.martonchik.o: ../../util/src/atmosphere.mfp.rayleigh.martonchik.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.rayleigh.martonchik.c


../../util/obj/atmosphere.mfp.spectralmodel.kato.o: ../../util/src/atmosphere.mfp.spectralmodel.kato.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.spectralmodel.kato.c


../../util/obj/atmosphere.mfp.spectralmodel.mlawer_lw.o: ../../util/src/atmosphere.mfp.spectralmodel.mlawer_lw.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.spectralmodel.mlawer_lw.c


../../util/obj/rrtm_driver_setcoef_taumol.o: ../../util/src/rrtm_driver_setcoef_taumol.f
	$(F77) $(FFLAGS) ../../util/src/rrtm_driver_setcoef_taumol.f


../../util/obj/setcoef.o: ../../util/src/setcoef.f
	$(F77) $(FFLAGS) ../../util/src/setcoef.f


../../util/obj/taumol.o: ../../util/src/taumol.f
	$(F77) $(FFLAGS) ../../util/src/taumol.f


../../util/obj/k_g.o: ../../util/src/k_g.f
	$(F77) $(FFLAGS) ../../util/src/k_g.f


../../util/obj/atmosphere.mfp.spectralmodel.pollack.o: ../../util/src/atmosphere.mfp.spectralmodel.pollack.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.spectralmodel.pollack.c


../../util/obj/atmosphere.mfp.template.o: ../../util/src/atmosphere.mfp.template.c
	 $(SUNCC) ../../util/src/atmosphere.mfp.template.c


../../util/obj/atmosphere.read.mcclatchey.o: ../../util/src/atmosphere.read.mcclatchey.c
	 $(SUNCC) ../../util/src/atmosphere.read.mcclatchey.c


../../util/obj/cntnm.xsect.o: ../../util/src/cntnm.xsect.c
	 $(SUNCC) ../../util/src/cntnm.xsect.c


../../util/obj/cntnmmodel.check.o: ../../util/src/cntnmmodel.check.c
	 $(SUNCC) ../../util/src/cntnmmodel.check.c


../../util/obj/cntnmmodel.read.o: ../../util/src/cntnmmodel.read.c
	 $(SUNCC) ../../util/src/cntnmmodel.read.c


../../util/obj/numrec.nrutil.o: ../../util/src/numrec.nrutil.c
	 $(SUNCC) ../../util/src/numrec.nrutil.c


../../util/obj/numrec.double.piksrt.o: ../../util/src/numrec.double.piksrt.c
	 $(SUNCC) ../../util/src/numrec.double.piksrt.c


../../util/obj/parsers.parameters.o: ../../util/src/parsers.parameters.c
	 $(SUNCC) ../../util/src/parsers.parameters.c


../../util/obj/phasefcns.henyeygreenstein.o: ../../util/src/phasefcns.henyeygreenstein.c
	 $(SUNCC) ../../util/src/phasefcns.henyeygreenstein.c


../../util/obj/phasefcns.read.ascii.o: ../../util/src/phasefcns.read.ascii.c
	 $(SUNCC) ../../util/src/phasefcns.read.ascii.c


../../util/obj/phasefcns.read.binary.o: ../../util/src/phasefcns.read.binary.c
	 $(SUNCC) ../../util/src/phasefcns.read.binary.c


../../util/obj/raprad.alloc.o: ../../util/src/raprad.alloc.c
	 $(SUNCC) ../../util/src/raprad.alloc.c


../../util/obj/raprad.atmosphere.layers.o: ../../util/src/raprad.atmosphere.layers.c
	 $(SUNCC) ../../util/src/raprad.atmosphere.layers.c


../../util/obj/raprad.atmosphere.layers.tau.o: ../../util/src/raprad.atmosphere.layers.tau.c
	 $(SUNCC) ../../util/src/raprad.atmosphere.layers.tau.c


../../util/obj/raprad.check.brdf.selection.o: ../../util/src/raprad.check.brdf.selection.c
	 $(SUNCC) ../../util/src/raprad.check.brdf.selection.c


../../util/obj/raprad.check.atmosphere_layers.o: ../../util/src/raprad.check.atmosphere_layers.c
	 $(SUNCC) ../../util/src/raprad.check.atmosphere_layers.c


../../util/obj/raprad.check.atmosphere_layers.lw.o: ../../util/src/raprad.check.atmosphere_layers.lw.c
	 $(SUNCC) ../../util/src/raprad.check.atmosphere_layers.lw.c


../../util/obj/raprad.check.constituents.o: ../../util/src/raprad.check.constituents.c
	 $(SUNCC) ../../util/src/raprad.check.constituents.c


../../util/obj/raprad.check.geometry.sunzenith.o: ../../util/src/raprad.check.geometry.sunzenith.c
	 $(SUNCC) ../../util/src/raprad.check.geometry.sunzenith.c


../../util/obj/raprad.check.photon.partition.o: ../../util/src/raprad.check.photon.partition.c
	 $(SUNCC) ../../util/src/raprad.check.photon.partition.c


../../util/obj/raprad.check.photon.space.final.o: ../../util/src/raprad.check.photon.space.final.c
	 $(SUNCC) ../../util/src/raprad.check.photon.space.final.c


../../util/obj/raprad.photon.space.final.o: ../../util/src/raprad.photon.space.final.c
	 $(SUNCC) ../../util/src/raprad.photon.space.final.c


../../util/obj/raprad.read.brdf.selection.o: ../../util/src/raprad.read.brdf.selection.c
	 $(SUNCC) ../../util/src/raprad.read.brdf.selection.c


../../util/obj/raprad.read.configuration.o: ../../util/src/raprad.read.configuration.c
	 $(SUNCC) ../../util/src/raprad.read.configuration.c


../../util/obj/raprad.read.constituent.phasefcn.o: ../../util/src/raprad.read.constituent.phasefcn.c
	 $(SUNCC) ../../util/src/raprad.read.constituent.phasefcn.c


../../util/obj/raprad.read.constituents.o: ../../util/src/raprad.read.constituents.c
	 $(SUNCC) ../../util/src/raprad.read.constituents.c


../../util/obj/raprad.read.geometry.sunzenith.o: ../../util/src/raprad.read.geometry.sunzenith.c
	 $(SUNCC) ../../util/src/raprad.read.geometry.sunzenith.c


../../util/obj/raprad.read.photon.partition.o: ../../util/src/raprad.read.photon.partition.c
	 $(SUNCC) ../../util/src/raprad.read.photon.partition.c


../../util/obj/raprad.read.photon.space.o: ../../util/src/raprad.read.photon.space.c
	 $(SUNCC) ../../util/src/raprad.read.photon.space.c


../../util/obj/rayleigh.bext.o: ../../util/src/rayleigh.bext.c
	 $(SUNCC) ../../util/src/rayleigh.bext.c


../../util/obj/rayleigh.crosssection.o: ../../util/src/rayleigh.crosssection.c
	 $(SUNCC) ../../util/src/rayleigh.crosssection.c


../../util/obj/rayleigh.indexrefraction.air.o: ../../util/src/rayleigh.indexrefraction.air.c
	 $(SUNCC) ../../util/src/rayleigh.indexrefraction.air.c


../../util/obj/rayleigh.meanfreepath.o: ../../util/src/rayleigh.meanfreepath.c
	 $(SUNCC) ../../util/src/rayleigh.meanfreepath.c


../../util/obj/rayleigh.opticaldepth.o: ../../util/src/rayleigh.opticaldepth.c
	 $(SUNCC) ../../util/src/rayleigh.opticaldepth.c


../../util/obj/spectralmodel.check.kato.o: ../../util/src/spectralmodel.check.kato.c
	 $(SUNCC) ../../util/src/spectralmodel.check.kato.c


../../util/obj/spectralmodel.check.mlawer_lw.o: ../../util/src/spectralmodel.check.mlawer_lw.c
	 $(SUNCC) ../../util/src/spectralmodel.check.mlawer_lw.c


../../util/obj/spectralmodel.check.pollack.o: ../../util/src/spectralmodel.check.pollack.c
	 $(SUNCC) ../../util/src/spectralmodel.check.pollack.c


../../util/obj/spectralmodel.read.kato.o: ../../util/src/spectralmodel.read.kato.c
	 $(SUNCC) ../../util/src/spectralmodel.read.kato.c


../../util/obj/spectralmodel.read.mlawer_lw.o: ../../util/src/spectralmodel.read.mlawer_lw.c
	 $(SUNCC) ../../util/src/spectralmodel.read.mlawer_lw.c


../../util/obj/spectralmodel.read.pollack.o: ../../util/src/spectralmodel.read.pollack.c
	 $(SUNCC) ../../util/src/spectralmodel.read.pollack.c


../../util/obj/thermodynamics.clausiusclapyron.o: ../../util/src/thermodynamics.clausiusclapyron.c
	 $(SUNCC) ../../util/src/thermodynamics.clausiusclapyron.c


../../util/obj/thermodynamics.units.o: ../../util/src/thermodynamics.units.c
	 $(SUNCC) ../../util/src/thermodynamics.units.c


raprad.2stream.lw.main.o: ../src/raprad.2stream.lw.main.c
	 $(SUNCC) ../src/raprad.2stream.lw.main.c


raprad.2stream.lw.add.o: ../src/raprad.2stream.lw.add.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.lw.add.f


raprad.2stream.lw.driver.o: ../src/raprad.2stream.lw.driver.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.lw.driver.f


raprad.2stream.lw.oppr1.o: ../src/raprad.2stream.lw.oppr1.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.lw.oppr1.f


raprad.2stream.lw.twostr.o: ../src/raprad.2stream.lw.twostr.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.lw.twostr.f


raprad.2stream.lw.newflux1.o: ../src/raprad.2stream.lw.newflux1.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.lw.newflux1.f


raprad.2stream.sw.main.o: ../src/raprad.2stream.sw.main.c
	 $(SUNCC) ../src/raprad.2stream.sw.main.c


raprad.2stream.sw.add.o: ../src/raprad.2stream.sw.add.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.sw.add.f


raprad.2stream.sw.driver.o: ../src/raprad.2stream.sw.driver.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.sw.driver.f


raprad.2stream.sw.twostr.o: ../src/raprad.2stream.sw.twostr.f
	$(F77) $(FFLAGS) ../src/raprad.2stream.sw.twostr.f


raprad.2stream.lw_mix.main.o: ../src_lw_mix/raprad.2stream.lw_mix.main.c
	 $(SUNCC) ../src_lw_mix/raprad.2stream.lw_mix.main.c


raprad.2stream.lw_mix.o: ../src_lw_mix/raprad.2stream.lw_mix.f
	$(F77) $(FFLAGS) ../src_lw_mix/raprad.2stream.lw_mix.f



clean:
	 rm $(OBJECTS) $(2STREAMLWOBJECTS) $(2STREAMSWOBJECTS)


#
#          END OF MAKE
#

