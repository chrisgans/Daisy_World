#!/usr/bin/env python3
# -*- coding: utf-8 -*-

if __name__ == "__main__":
    print("In __main__")
    import numpy as np
    import matplotlib.pyplot as plt

    import model_functions as mf

    ### 4)
    luminosities = np.arange(0.5, 1.6, 0.002) # Stelar luminosities
    alphaw_out = np.ones(len(luminosities)) * np.nan # Output variable for white
    alphab_out = np.ones(len(luminosities)) * np.nan # Output variable for black
    temp_out = np.ones(len(luminosities)) * np.nan   # Output variable for temp

    gamma_arr = [0.1, 0.2, 0.3, 0.5, 0.7]
    fig, axs = plt.subplots(len(gamma_arr), 1, figsize=(10, 15))

    for k, gamma in enumerate(gamma_arr):
      # Define constants and variables
      alphaw = mf.alphaw   # Cover fraction of white daisies
      alphab = mf.alphab   # Cover fraction of black daisies

      ### custom constants
      alphaw_min = mf.alphaw
      alphab_min = mf.alphab

      # Main loop for changing luminosity
      for i,L in enumerate(luminosities):
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
              planet_albedo = mf.albedo(alphaw, alphab, alphag, mf.aw, mf.ab, mf.ag)
              # Planetary temperature
              T = mf.planetary_temp(mf.S, planet_albedo, L=L)
              # Local temperature
              Tw = mf.local_temp(planet_albedo, mf.aw, T)
              Tb = mf.local_temp(planet_albedo, mf.ab, T)
              # Birth rate
              betaw = mf.beta(Tw)
              betab = mf.beta(Tb)
              # Change in daisies
              dawdt = mf.daisy_replicator(alphaw, alphag, betaw, gamma)
              dabdt = mf.daisy_replicator(alphab, alphag, betab, gamma)
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
      axs[k].plot(luminosities, alphaw_out*100, 'b', label='White')
      axs[k].plot(luminosities, alphab_out*100, 'k', label='Black')
      axs[k].set_ylabel('gamma = %s' % gamma)

      if (k == 0):
        axs[0].legend(loc='upper right')

    plt.tight_layout()
    plt.savefig('./outputs/exc_4.png', dpi=520)
