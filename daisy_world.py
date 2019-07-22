#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 19:19:16 2019

@author: tammas

This is an implementation of daisy world: a theoretical
model of the Gaia hypothesis. This theory was initially developed by
Lovelock (1983) to demonstrate the plausibility of living things
interacting with, and regulating, their environment.

    Wood, A. J., G. J. Ackland, J. G. Dyke, H. T. P. Williams, and T. M.
        Lenton, 2015: Daisywolrd: a Review. Rev. Geophys., 46, RG1001,
        https://doi.org/10.1029/2006RG000217.
"""


import numpy as np
import matplotlib.pyplot as plt


def albedo(alphaw, alphab, alphag, aw, ab, ag):
    """Caluclate the average planetary albedo.
    alphaw*aw + alphab*ab + alphag*ag

    Arguments
    ---------
    alphaw : float
        Surface fraction of white daisies.
    alphab : float
        Surface fraction of black daisies.
    alphag : float
        Surface fraction of bare ground.
    aw = 0.7 : float
        Albedo of white daisies.
    ab = 0.8 : float
        Albedo of black dasies.
    ag = 0.5 : float
        Albedo of bare ground
    """
    return alphaw*aw + alphab*ab + alphag*ag


def daisy_replicator(alpha, alphag, beta, gamma):
    """Calculate the rate of change of a replicator (daisies). The formula is:
    a*[ag*b-g]
    which describes a logistic growth and constant death within a limited
    resource system over time t.

    Arguments
    ---------
    alpha :  float
        alpha -- the fraction of surface covered by daisies.
    alphag : float
        alpha_g -- The fraction of available bare ground.
    beta :  float
        beta -- the birth rate.
    gamma : float
        gamma -- the death rate.
    """
    return alpha*(alphag*beta-gamma)


def beta(temp, optimum=295.65, tmin=278.15, tmax=313.15, k=0.003265):
    """Calculate the environment parameter beta. This is like a birth rate
    that reaches a maximum at an optimal temperature, and a range of birth
    rates is specified by a parabolic width parameter of k=0.003265.
    1 - k*(temp-optimum)**2

    Arguments
    ---------
    temp : float
        The environment temperature experienced by the daisies

    Keyword arguments
    -----------------
    optimum = 295.65 : float
        The optimum temperature for daisy birthrate.
    tmin = 278.15 : float
        Maximum temperature for birthrate.
    tmax = 313.15 : float
        Minimum temperature for birthrate.
    k = 0.003265 : float
        Quadratic width parameter.
    """
    if (temp>tmin)&(temp<tmax):
        return 1 - k*(temp-optimum)**2
    else:
        return 0


def euler(initial, tendency, h=1):
    """Integrate forward in time using Euler's method of numerical integration.
    initial + h*tendency

    Arguments
    ---------
    initial :  float
        The initial state.
    tendency : float
        The rate of change in the initial state.

    Keyword arguments
    -----------------
    h = 1 : float
        The timestep duration.
    """
    return initial + h*tendency


def local_temp(A, albedo, T, q=30):
    """Calculate local temperature experienced by a particular daisy type or
    ground cover. This is a simplified version of the original.
    q*(A-albedo)+T

    Arguments
    ---------
    A : float
        Planetary albedo
    alpha : float
        Albedo of daisy type
    T : float
        Planetary temperature.

    Keyword Arguments
    -----------------
    q = 30 : float
        Heat transfer coefficient
    """
    return q*(A-albedo)+T


def planetary_temp(S, A, L=1.0):
    """Calculate the planetary temperature.
    SL(1-A) = sT**4

    Arguments
    ---------
    S : float
        Incident solar energy.
    A : float
        Planetary albedo.

    Keyword Arguments
    -----------------
    L = 1.0 : float
        Normalised stellar luminosity.
    """
    sigma = 5.67032e-8 # Stephan-Bolzmann constant.
    return ((S*L*(1-A))/sigma)**(1/4.)


### 1)
fig, axs = plt.subplots(2, 3)

for e, (w, b) in enumerate(zip([0.01, 0.00, 0.01], [0.00, 0.01, 0.01])):
  # Define constants and variables
  alphaw = w   # Cover fraction of white daisies
  alphab = b   # Cover fraction of black daisies
  p = 1           # The fraction of habitable surface
  alphag = p-alphaw-alphab # Cover fraction of bare ground
  aw = 0.75       # Albedo of white daisies
  ab = 0.25       # Albedo of black daisies
  ag = 0.5        # Albedo of bare ground
  gamma = 0.3     # The death rate
  S = 1000        # Solar constant (W/m^2)
  maxn = 1000     # Maximum number of iterations
  tol = 0.000001  # Tollerance of solution
  luminosities = np.arange(0.5,1.6, 0.002) # Stelar luminosities
  alphaw_out = np.ones(len(luminosities)) * np.nan # Output variable for white
  alphab_out = np.ones(len(luminosities)) * np.nan # Output variable for black
  temp_out = np.ones(len(luminosities)) * np.nan   # Output variable for temp


  ### custom constants
  alphaw_min = w
  alphab_min = b

  # Main loop for changing luminosity
  for i,L in enumerate(luminosities):
      # Set a minimum for cover fractions
      alphaw = alphaw_min if alphaw < alphaw_min else alphaw
      alphab = alphab_min if alphab < alphab_min else alphab
      alphag = p-alphaw-alphab
      # Reset counters
      n = 0
      changew, changeb = 1,1
      # Run loop for daisy earth.
      while (n<maxn) and (changew>tol) and (changeb>tol):
          # Store the initial cover fractions
          sw,sb = alphaw, alphab
          # Planetary albedo
          planet_albedo = albedo(alphaw,alphab,alphag,aw,ab,ag)
          # Planetary temperature
          T = planetary_temp(S,planet_albedo, L=L)
          # Local temperature
          Tw = local_temp(planet_albedo,aw,T)
          Tb = local_temp(planet_albedo,ab,T)
          # Birth rate
          betaw = beta(Tw)
          betab = beta(Tb)
          # Change in daisies
          dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
          dabdt = daisy_replicator(alphab, alphag, betab, gamma)
          # Integrate
          alphaw = euler(alphaw, dawdt)
          alphab = euler(alphab, dabdt)
          alphag = p-alphaw-alphab
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
  axs[1, e].plot(luminosities, temp_out-273.15,'r')
  if e == 0:
    axs[1, e].set_ylabel('Temperature (°C)')
  if e == 1:
    axs[1, e].set_xlabel('Luminosity')
    axs[1, e].set_title('Planetary temperature')
  plt.show()




### 2)
aw_arr = [0.00, 0.25, 0.75, 0.9, 1]
ab_arr = [0.00, 0.25, 0.75, 0.9, 1]
fig, axs = plt.subplots(len(aw_arr), len(ab_arr))

for k, aw in enumerate(aw_arr):
  for j, ab in enumerate(ab_arr):
    # Define constants and variables
    alphaw = 0.01   # Cover fraction of white daisies
    alphab = 0.01   # Cover fraction of black daisies
    p = 1           # The fraction of habitable surface
    alphag = p-alphaw-alphab # Cover fraction of bare ground
    #  aw = 0.75       # Albedo of white daisies
    #  ab = 0.25       # Albedo of black daisies
    ag = 0.5        # Albedo of bare ground
    gamma = 0.3     # The death rate
    S = 1000        # Solar constant (W/m^2)
    maxn = 1000     # Maximum number of iterations
    tol = 0.000001  # Tollerance of solution
    luminosities = np.arange(0.5,1.6, 0.002) # Stelar luminosities
    alphaw_out = np.ones(len(luminosities)) * np.nan # Output variable for white
    alphab_out = np.ones(len(luminosities)) * np.nan # Output variable for black
    temp_out = np.ones(len(luminosities)) * np.nan   # Output variable for temp


    ### custom constants
    alphaw_min = 0.01
    alphab_min = 0.01

    # Main loop for changing luminosity
    for i,L in enumerate(luminosities):
        # Set a minimum for cover fractions
        alphaw = alphaw_min if alphaw < alphaw_min else alphaw
        alphab = alphab_min if alphab < alphab_min else alphab
        alphag = p-alphaw-alphab
        # Reset counters
        n = 0
        changew, changeb = 1,1
        # Run loop for daisy earth.
        while (n<maxn) and (changew>tol) and (changeb>tol):
            # Store the initial cover fractions
            sw,sb = alphaw, alphab
            # Planetary albedo
            planet_albedo = albedo(alphaw, alphab, alphag, aw, ab, ag)
            # Planetary temperature
            T = planetary_temp(S,planet_albedo, L=L)
            # Local temperature
            Tw = local_temp(planet_albedo,aw,T)
            Tb = local_temp(planet_albedo,ab,T)
            # Birth rate
            betaw = beta(Tw)
            betab = beta(Tb)
            # Change in daisies
            dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
            dabdt = daisy_replicator(alphab, alphag, betab, gamma)
            # Integrate
            alphaw = euler(alphaw, dawdt)
            alphab = euler(alphab, dabdt)
            alphag = p-alphaw-alphab
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

    plt.show()


### 3)
w_opt_arr = [278.15, 285.0, 295.65, 305.0, 313.15]
b_opt_arr = [278.15, 285.0, 295.65, 305.0, 313.15]
fig, axs = plt.subplots(len(w_opt_arr), len(b_opt_arr))

for k, w_opt in enumerate(w_opt_arr):
  for j, b_opt in enumerate(b_opt_arr):
    # Define constants and variables
    alphaw = 0.01   # Cover fraction of white daisies
    alphab = 0.01   # Cover fraction of black daisies
    p = 1           # The fraction of habitable surface
    alphag = p-alphaw-alphab # Cover fraction of bare ground
    aw = 0.75       # Albedo of white daisies
    ab = 0.25       # Albedo of black daisies
    ag = 0.5        # Albedo of bare ground
    gamma = 0.3     # The death rate
    S = 1000        # Solar constant (W/m^2)
    maxn = 1000     # Maximum number of iterations
    tol = 0.000001  # Tollerance of solution
    luminosities = np.arange(0.5, 1.6, 0.002) # Stelar luminosities

    alphaw_out = np.ones(len(luminosities)) * np.nan # Output variable for white
    alphab_out = np.ones(len(luminosities)) * np.nan # Output variable for black
    temp_out = np.ones(len(luminosities)) * np.nan   # Output variable for temp


    ### custom constants
    alphaw_min = 0.01
    alphab_min = 0.01

    # Main loop for changing luminosity
    for i,L in enumerate(luminosities):
        # Set a minimum for cover fractions
        alphaw = alphaw_min if alphaw < alphaw_min else alphaw
        alphab = alphab_min if alphab < alphab_min else alphab
        alphag = p-alphaw-alphab
        # Reset counters
        n = 0
        changew, changeb = 1,1
        # Run loop for daisy earth.
        while (n<maxn) and (changew>tol) and (changeb>tol):
            # Store the initial cover fractions
            sw,sb = alphaw, alphab
            # Planetary albedo
            planet_albedo = albedo(alphaw, alphab, alphag, aw, ab, ag)
            # Planetary temperature
            T = planetary_temp(S,planet_albedo, L=L)
            # Local temperature
            Tw = local_temp(planet_albedo,aw,T)
            Tb = local_temp(planet_albedo,ab,T)
            # Birth rate
            betaw = beta(Tw, optimum=w_opt)
            betab = beta(Tb, optimum=b_opt)
            # Change in daisies
            dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
            dabdt = daisy_replicator(alphab, alphag, betab, gamma)
            # Integrate
            alphaw = euler(alphaw, dawdt)
            alphab = euler(alphab, dabdt)
            alphag = p-alphaw-alphab
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
      axs[k, j].set_ylabel('w_opt = %s' % w_opt)
    if k == len(b_opt_arr) - 1:
      axs[k, j].set_xlabel('b_opt = %s' % b_opt)

    plt.show()


### 4)
gamma_arr = [0.1, 0.2, 0.3, 0.5, 0.7]
fig, axs = plt.subplots(len(gamma_arr), 1)

for k, gamma in enumerate(gamma_arr):
  # Define constants and variables
  alphaw = 0.01   # Cover fraction of white daisies
  alphab = 0.01   # Cover fraction of black daisies
  p = 1           # The fraction of habitable surface
  alphag = p-alphaw-alphab # Cover fraction of bare ground
  aw = 0.75       # Albedo of white daisies
  ab = 0.25       # Albedo of black daisies
  ag = 0.5        # Albedo of bare ground
  ### gamma = 0.3     # The death rate
  S = 1000        # Solar constant (W/m^2)
  maxn = 1000     # Maximum number of iterations
  tol = 0.000001  # Tollerance of solution
  luminosities = np.arange(0.5, 1.6, 0.002) # Stelar luminosities

  alphaw_out = np.ones(len(luminosities)) * np.nan # Output variable for white
  alphab_out = np.ones(len(luminosities)) * np.nan # Output variable for black
  temp_out = np.ones(len(luminosities)) * np.nan   # Output variable for temp


  ### custom constants
  alphaw_min = 0.01
  alphab_min = 0.01

  # Main loop for changing luminosity
  for i,L in enumerate(luminosities):
      # Set a minimum for cover fractions
      alphaw = alphaw_min if alphaw < alphaw_min else alphaw
      alphab = alphab_min if alphab < alphab_min else alphab
      alphag = p-alphaw-alphab
      # Reset counters
      n = 0
      changew, changeb = 1,1
      # Run loop for daisy earth.
      while (n<maxn) and (changew>tol) and (changeb>tol):
          # Store the initial cover fractions
          sw,sb = alphaw, alphab
          # Planetary albedo
          planet_albedo = albedo(alphaw, alphab, alphag, aw, ab, ag)
          # Planetary temperature
          T = planetary_temp(S,planet_albedo, L=L)
          # Local temperature
          Tw = local_temp(planet_albedo,aw,T)
          Tb = local_temp(planet_albedo,ab,T)
          # Birth rate
          betaw = beta(Tw)
          betab = beta(Tb)
          # Change in daisies
          dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
          dabdt = daisy_replicator(alphab, alphag, betab, gamma)
          # Integrate
          alphaw = euler(alphaw, dawdt)
          alphab = euler(alphab, dabdt)
          alphag = p-alphaw-alphab
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
  plt.show()





### 5)
p = 1           # The fraction of habitable surface
alphag = p - alphaw - alphab # Cover fraction of bare ground
aw = 0.75       # Albedo of white daisies
ab = 0.25       # Albedo of black daisies
ag = 0.5        # Albedo of bare ground
gamma = 0.3     # The death rate
S = 1000        # Solar constant (W/m^2)
maxn = 1000     # Maximum number of iterations
tol = 0.000001  # Tollerance of solution


### 'custom constants'
L = 1
alphaw = 0.2
alphab = 0.5
alphag = p-alphaw-alphab
sw, sb = 0, 0
n = 0

alphaw_out = [alphaw]
alphab_out = [alphab]
planet_temp_out = []

while (n < maxn) & ((abs(alphaw - sw) > 0.001) | (abs(alphab - sb) > 0.001)):
  # Store the initial cover fractions
  sw, sb = alphaw, alphab
  # Planetary albedo
  planet_albedo = albedo(alphaw, alphab, alphag, aw, ab, ag)
  # Planetary temperature
  T = planetary_temp(S,planet_albedo, L=L)
  planet_temp_out.append(T)
  # Local temperature
  Tw = local_temp(planet_albedo, aw, T)
  Tb = local_temp(planet_albedo, ab, T)
  # Birth rate
  betaw = beta(Tw)
  betab = beta(Tb)
  # Change in daisies
  dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
  dabdt = daisy_replicator(alphab, alphag, betab, gamma)
  # Integrate
  alphaw = euler(alphaw, dawdt)
  alphab = euler(alphab, dabdt)

  alphaw_out.append(alphaw)
  alphab_out.append(alphab)

  alphag = p - alphaw - alphab
  n += 1


fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(alphaw_out, "b-", label="Alpha_w")
ax1.plot(alphab_out, "k-", label="Alpha_b")
ax1.set_ylabel("Cover Fraction in %")
ax1.legend()

ax2.plot(np.array(planet_temp_out) - 273.15, "r-")
ax2.set_ylabel("Temperature in C°")
ax2.set_xlabel("Iteration")

plt.show()

