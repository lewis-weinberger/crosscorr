#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to plot power spectrum output of crosscorr
"""

import sys
import numpy as np
import matplotlib.pyplot as plt


if (len(sys.argv) != 2):
    print("Usage: python3 plot.py [OUTPUT FILE]")
    sys.exit()

filename = sys.argv[1].strip()

k, powerspec, deltasqk, _ = np.loadtxt(filename, unpack=True)
        
print('{} read in'.format(filename))

################################################################################

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1,1,1)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(r'$k$ [1/units]')
ax.set_ylabel(r'$P(k)$')

x = k
y = np.absolute(powerspec)

ax.plot(x, y, ls='-', c='r', lw=3)
ax.tick_params(axis='both', which='major', pad=15)

fig.tight_layout()
plt.savefig("crosscorr_1.png", bbox_inches="tight")

################################################################################


fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1,1,1)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel(r'$k$ [1/units]')
ax.set_ylabel(r'$|\Delta^2|(k)$')

x = k
y = np.absolute(deltasqk)

ax.plot(x, y, ls='-', c='r', lw=3)

ax.tick_params(axis='both', which='major', pad=15)

fig.tight_layout()
plt.savefig("crosscorr_2.png", bbox_inches="tight")
