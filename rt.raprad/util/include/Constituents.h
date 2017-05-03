typedef struct CONSTITUENTS {

	char
		name[64],
		phasefcn_file[128],
		phasefcn_format[8],
		mfp_fcn_name[64],
		es[8];

	double
		base_height,
		top_height,
		scale_height;
	int
		base_height_index,
		top_height_index;

	int
		p_which,
		p_explicit_number,
		p_legendre_number;

	double
		cdensity,
		tau,
		*tau_user,
		tau_sca,
		*tau_sca_user,
		tau_abs,
		*tau_abs_user,
		w0,
		*w0_user,
		sig_ext,
		*sig_ext_user,
		sig_sca,
		*sig_sca_user,
		sig_abs,
		*sig_abs_user,
		k_ext,
		k_sca,
		k_abs;

	double
		p_dfcn_strength_o,
		p_dfcn_strength_f,
		p_normalization,
		p_g,
		*p_g_user,
		p_g_file,
		p_g_dfcn,
		*p_legendre_coef,
		*pmu,
		*p,
		*pdfcno,
		*pdfcnf,
		**pqn,
		**pqp,
		*psqn,
		*psqp,
		*psvn,
		*psvp;

	float
		*costheta,
		*pscat,
		*derivpen;


} Constituents;
