
#define TINY 1.E-10

typedef struct RADIATIVE_TRANSFER_WORK {

	int
		*indx_zxa;

	float

		**rint,
		**deriv,

		*rint_z,
		*deriv_z;

	double

		*col_zxa,

		taua,
		taura,
		taub,
		taurb,

		**r,
		**t,
		**rntot,
		**tptot,

		**tr,
		**rl,
		**rr,

		**tmp,

		**rnu,
		**rpu,
		**tnu,
		**tpu,

		**rn,
		**rp,
		**tn,
		**tp,

		*Y,
		*Ytmp,
		*Z,
		*Ztmp,

		*bbradlayern,
		*bbradlayerp,
		*bbrad,
		*bbradn,
		*bbradp,

		*br,
		*bt,

		**layers_ss_pqn,
		**layers_ss_pqp,
		*layers_ss_psqn,
		*layers_ss_psqp,
		*layers_ss_psvn,
		*layers_ss_psvp,

		***raddown,
		**radup;

} RtWork;

typedef struct RADIATIVE_TRANSFER {

	int
		layers_number,
		rt_number,
		layer_doubling_number,
		*layers_scatterflag;

	float
		layer_doubling_tau,
		**layers_kext_c,
		**layers_ksca_c,
		**layers_kabs_c,
		**layers_tau_c,
		**layers_tau_c_dfcn_o,
		**layers_tau_c_dfcn_f,

		*layers_tau,
		*layers_tau_dfcn_o,
		*layers_tau_dfcn_f,

		*layers_taucum,
		*layers_taucum_dfcn_o,
		*layers_taucum_dfcn_f,

		*layers_kext,
		*layers_kext_dfcn_o,
		*layers_kext_dfcn_f,

		*layers_w0,
		*layers_w0_dfcn_o,
		*layers_w0_dfcn_f,

		***layers_pqn,
		***layers_pqp,

                *layers_legendre_coef,
                *layers_g0;

	double
		*layers_bb,

		*plankbnd,
		*planklayd,
		*planklayu,
		*planklevd,
		*planklevu,

		****tn,
		****tp,
		****rn,
		****rp,

		 *fluxp,
		 *fluxn,

		**radsqp,
		**radsqn,
		**radsvp,
		**radsvn,

		  radsundir,
		  radsundirdelta,

		 *radaprsqp,
		 *radaprsqn,
		 *radaprsvp,
		 *radaprsvn,

		 *radaprdirsqp,
		 *radaprdirsvp,
		 *radaprdirdeltasqp,
		 *radaprdirdeltasvp,

		**radaprssqp,
		**radaprssqn,

		 *radaprsssqp,
		 *radaprsssqn,
		 *radaprsssvp,
                 *radaprsssvn,

		 *radaprmltsqp,
		 *radaprmltsqn,
		 *radaprmltsvp,
		 *radaprmltsvn,

		 **radsrfsqp,
		 **radsrfsqn,
		 **radsrfsvp,
		 **radsrfsvn,

		***radsrfdirsqn,
		***radsrfdirsvn,
		***radsrfdirdeltasqn,
		***radsrfdirdeltasvn,

		***radsrfmltsqp,
		***radsrfmltsqn,
		***radsrfmltsvp,
		***radsrfmltsvn;

} Rt;
