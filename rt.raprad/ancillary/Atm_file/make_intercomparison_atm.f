      program make_intercomparison_atm

% This file generates atm file from intercomparison profile.
      parameter(n_layer = 100)

      real level_no(n_layer), level_km(n_layer), press_mb(n_layer)
      real temp_k(n_layer), ozone_g_g(n_layer)
      real layer_no(n_layer), mid_layer_km(n_layer)
      real density_kg_m3(n_layer), cf_liq(n_layer)
      real cf_ice(n_layer), cf_both(n_layer)      
      real acf_dn_liq(n_layer), acf_dn_ice(n_layer)
      real acf_dn_both(n_layer), acf_up_liq(n_layer), 
      real acf_up_ice(n_layer), acf_up_both(n_layer)    
      real mr_liq(n_layer), mr_ice(n_layer), mr_both(n_layer)
      real std_mr_liq(n_layer), std_mr_ice(n_layer)
      real std_mr_both(n_layer), log_mr_liq(n_layer)    
      real log_mr_ice(n_layer), log_mr_both(n_layer)
      real wv_mr_clr(n_layer)     
      real wv_mr_cld(n_layer), wv_mr_all(n_layer)



      open(unit=11,file='intercomparison.original.atm.dat'
     &     status='old')

      do 1000 i = 1, 65
        read(11,*) level_no(i), level_km(i), press_mb(i),
     &             temp_k(i), ozone_g_g(i)
        read(11,100)
        read(11,*) layer_no(i), mid_layer_km(i), density_kg_m3(i), 
     &             cf_liq(i), cf_ice(i)
        read(11,*) cf_both(i), acf_dn_liq(i), acf_dn_ice(i),
     &             acf_dn_both(i), acf_up_liq(i)
        read(11,*) acf_up_ice(i), acf_up_both(i), mr_liq(i), 
     &             mr_ice(i), mr_both(i)
        read(11,*) std_mr_liq(i), std_mr_ice(i), std_mr_both(i),   
     &             log_mr_liq(i), log_mr_ice(i)
        read(11,*) log_mr_both(i), wv_mr_clr(i),
     &             wv_mr_cld(i), wv_mr_all(i)
 1000 continue 

      i = 66

      read(11,*) level_no(i), level_km(i), press_mb(i),
     &           temp_k(i), ozone_g_g(i)




  100 format(1x)
