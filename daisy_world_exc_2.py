#!/usr/bin/env python3
# -*- coding: utf-8 -*-

if __name__ == "__main__":
    print("In __main__")
    import numpy as np
    import matplotlib.pyplot as plt

    import model_functions as mf

    ### 2)
    luminosities = np.arange(0.5,1.6, 0.002) # Stelar luminosities
    alphaw_out = np.ones(len(luminosities)) * np.nan # Output variable for white
    alphab_out = np.ones(len(luminosities)) * np.nan # Output variable for black
    temp_out = np.ones(len(luminosities)) * np.nan   # Output variable for temp

    aw_arr = [0.00, 0.25, 0.75, 0.9, 1]
    ab_arr = [0.00, 0.25, 0.75, 0.9, 1]
    fig, axs = plt.subplots(len(aw_arr), len(ab_arr), figsize=(15, 15))

    for k, aw in enumerate(aw_arr):
      for j, ab in enumerate(ab_arr):
        # Define constants and variables
        alphaw = mf.alphaw   # Cover fraction of white daisies
        alphab = mf.alphab   # Cover fraction of black daisies

        ### custom constants
        alphaw_min = mf.alphaw
        alphab_min = mf.alphab

        # Main loop for changing luminosity
        for i, L in enumerate(luminosities):
            # Set a minimum for cover fractions
            alphaw = alphaw_min if alphaw < alphaw_min else alphaw
            alphab = alphab_min if alphab < alphab_min else alphab
            alphag = mf.p - alphaw - alphab
            # Reset counters
            n = 0
            changew, changeb = 1,1
            # Run loop for daisy earth.
            while (n < mf.maxn) and (changew > mf.tol) and (changeb > mf.tol):
                # Store the initial cover fractions
                sw, sb = alphaw, alphab
                # Planetary albedo
                planet_albedo = mf.albedo(alphaw, alphab, alphag, aw, ab, mf.ag)
                # Planetary temperature
                T = mf.planetary_temp(mf.S, planet_albedo, L=L)
                # Local temperature
                Tw = mf.local_temp(planet_albedo, aw, T)
                Tb = mf.local_temp(planet_albedo, ab, T)
                # Birth rate
                betaw = mf.beta(Tw)
                betab = mf.beta(Tb)
                # Change in daisies
                dawdt = mf.daisy_replicator(alphaw, alphag, betaw, mf.gamma)
                dabdt = mf.daisy_replicator(alphab, alphag, betab, mf.gamma)
                # Integrate
                alphaw = mf.euler(alphaw, dawdt)
                alphab = mf.euler(alphab, dabdt)
                alphag = mf.p - alphaw - alphab
                n += 1
            # Store the output
            alphaw_out[i] = alphaw
            alphab_out[i] = alphab
            temp_out[i] = T

        # Plot the results
        # Cover fractions
        axs[k, j].plot(luminosities, alphaw_out*100, 'b', label='White')
        axs[k, j].plot(luminosities, alphab_out*100, 'k', label='Black')

        if (k == 0) & (j == 0):
          axs[0, 0].legend(loc='upper right')
        if j == 0:
          axs[k, j].set_ylabel('aw = %s' % aw)
        if k == len(aw_arr) - 1:
          axs[k, j].set_xlabel('ab = %s' % ab)

        plt.tight_layout()
        plt.savefig('./outputs/exc_2.png', dpi=340)
