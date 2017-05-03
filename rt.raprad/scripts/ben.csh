#!/bin/csh
# chmod +x run.csh to change the permission for the script file

set pathcode = "~/research/radiation/Raprad-9/rt.raprad/2stream/obj"
set pathflux = "~/research/radiation/Raprad-9/rt.raprad/rapradtally/obj"


cd $pathcode 
make clean
make || exit 1
raprad2streamlw ../../scenarios/bbflux/st_atm_lw/Rt1d.configuration.lw 0
raprad2streamsw ../../scenarios/bbflux/st_atm/Rt1d.configuration.all_H2O 0

cd $pathflux
make clean
make || exit 1
rapradbbfluxlw ../../scenarios/bbflux/st_atm_lw/Rt1d.configuration.lw 0
rapradbbfluxsw ../../scenarios/bbflux/st_atm/Rt1d.configuration.all_H2O 0

combine_sw_lw_flx
