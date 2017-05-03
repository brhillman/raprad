#!/bin/csh
# chmod +x run_all.csh to change the permission for the script file

set pathcode = "/home/disk/brume2/tra/Research/radiation/Tom/raprad/Raprad-8/rt.raprad/2stream/obj"
set pathflux = "/home/disk/brume2/tra/Research/radiation/Tom/raprad/Raprad-8/rt.raprad/rapradtally/obj"


cd $pathcode 
make clean
make
raprad2streamlw ../../scenarios/bbflux/st_atm_lw/Rt1d.configuration.lw 0
#raprad2streamsw ../../scenarios/bbflux/st_atm/Rt1d.configuration.all_H2O 0

cd $pathflux
make clean
make
rapradbbfluxlw ../../scenarios/bbflux/st_atm_lw/Rt1d.configuration.lw 0
#rapradbbfluxsw ../../scenarios/bbflux/st_atm/Rt1d.configuration.all_H2O 0

#combine_sw_lw_flx