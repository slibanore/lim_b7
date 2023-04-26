import fisher_pmf_b7 as b7

def main():

    n_B = 2.9
    sigma_B_0 = 1.
    smooth_scale = 10.
    model_id = 'COS3_hz'
    derive_pars = ['n_B','sigma_B_0','astro_alpha','astro_beta','astro_sig_SFR', 'astro_sig_scatter','astro_Lcut']

    save_VID_flag = True
    save_fisher_flag = True
    import_allowed = True

    plot_VID_flag = True
    plot_sigma_flag = True
    
    prior_matrix = False
    
    b7.fisher_VID_pmf(n_B = n_B, 
                    sigma_B_0 = sigma_B_0, 
                    smooth_scale = smooth_scale,
                    model_id = model_id, \
                    derive_pars = derive_pars,\
                    save_VID_flag = save_VID_flag, \
                    save_fisher_flag = save_fisher_flag, \
                    plot_VID_flag = plot_VID_flag, \
                    plot_sigma_flag = plot_sigma_flag, \
                    import_allowed = import_allowed,
                    prior_matrix = prior_matrix)


if __name__ == "__main__":
    main()