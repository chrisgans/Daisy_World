#!/usr/bin/env python3
# -*- coding: utf-8 -*-

if __name__ == "__main__":
    print("In __main__")
    import numpy as np
    import matplotlib.pyplot as plt

    import model_functions as mf



    ### 5)
    ### 'custom constants'
    maxn = 1000     # Maximum number of iterations
    tol = 0.000001  # Tollerance of solution
    L = 1
    alphaw = 0.2
    alphab = 0.5
    alphag = mf.p - alphaw - alphab
    sw, sb = 0, 0
    n = 0

    alphaw_out = [alphaw]
    alphab_out = [alphab]
    planet_temp_out = []

    while (n < maxn) & ((abs(alphaw - sw) > 0.001) | (abs(alphab - sb) > 0.001)):
      # Store the initial cover fractions
      sw, sb = alphaw, alphab
      # Planetary albedo
      planet_albedo = mf.albedo(alphaw, alphab, alphag, mf.aw, mf.ab, mf.ag)
      # Planetary temperature
      T = mf.planetary_temp(mf.S, planet_albedo, L=L)
      planet_temp_out.append(T)
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

      alphaw_out.append(alphaw)
      alphab_out.append(alphab)

      alphag = mf.p - alphaw - alphab
      n += 1


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 15))
    ax1.plot(alphaw_out, "b-", label="Alpha_w")
    ax1.plot(alphab_out, "k-", label="Alpha_b")
    ax1.set_ylabel("Cover Fraction in %")
    ax1.legend()

    ax2.plot(np.array(planet_temp_out) - 273.15, "r-")
    ax2.set_ylabel(u"Temperature in C°")
    ax2.set_xlabel("Iteration")

    plt.tight_layout()
    plt.savefig('./outputs/exc_5.png', dpi=520)

