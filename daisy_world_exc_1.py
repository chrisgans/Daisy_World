#!/usr/bin/env python3
# -*- coding: utf-8 -*-

if __name__ == "__main__":
    print("In __main__")
    import numpy as np
    import matplotlib.pyplot as plt

    import model_functions as mf


    ### 1)
    luminosities = np.arange(0.5, 1.6, 0.002) # Stelar luminosities
    alphaw_out = np.ones(len(luminosities)) * np.nan # Output variable for white
    alphab_out = np.ones(len(luminosities)) * np.nan # Output variable for black
    temp_out = np.ones(len(luminosities)) * np.nan   # Output variable for temp

    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    for e, (w, b) in enumerate(zip([0.01, 0.00, 0.01], [0.00, 0.01, 0.01])):
      # Define variables
      alphaw = w   # Cover fraction of white daisies
      alphab = b   # Cover fraction of black daisies

      ### custom constants
      alphaw_min = w
      alphab_min = b

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
      axs[0, e].plot(luminosities, alphaw_out*100, 'b', label='White')
      axs[0, e].plot(luminosities, alphab_out*100, 'k', label='Black')
      if e == 0:
        axs[0, e].legend(loc='upper right')
        axs[0, e].set_ylabel('Surface cover %')

      if e == 1:
        axs[0, e].set_xlabel('Luminosity')
        axs[0, e].set_title('Cover fractions')

      # Planetary temperature
      axs[1, e].plot(luminosities, temp_out-273.15, 'r')
      if e == 0:
        axs[1, e].set_ylabel(u'Temperature (°C)')
      if e == 1:
        axs[1, e].set_xlabel('Luminosity')
        axs[1, e].set_title('Planetary temperature')

    plt.tight_layout()
    plt.savefig('./outputs/exc_1.png', dpi=520)
