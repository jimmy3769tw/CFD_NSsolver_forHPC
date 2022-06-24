#!/usr/bin/env python
from pylab import *
from scipy.io import mmread

A = mmread('This_Matrix.mtx')

fig, (ax1, ax2) = subplots(2, 1, sharex=True, figsize=(8,10), gridspec_kw=dict(height_ratios=[4,1]))
ax1.spy(A, marker='.', markersize=0.5, alpha=0.5)




axins = ax1.inset_axes([0.55, 0.55, 0.45, 0.45])
axins.spy(A, marker='o', markersize=3, alpha=0.5)
n = (A.shape[0] // 3 // 3) * 3
axins.set_xlim([n - 0.5, n + 44.5])
axins.set_ylim([n - 0.5, n + 44.5])
axins.invert_yaxis()
axins.set_xticklabels('')
axins.set_yticklabels('')
ax1.indicate_inset_zoom(axins)



ax2.spy(A, marker='.', markersize=0.5, alpha=0.5)

axins2 = ax2.inset_axes([0.55, 0.55, 0.45, 0.45])
axins2.spy(A, marker='o', markersize=3, alpha=0.5)
axins2.set_xlim([n - 0.5, n + 44.5])
axins2.set_ylim([n - 0.5, n + 44.5])
axins2.invert_yaxis()
axins2.set_xticklabels('')
axins2.set_yticklabels('')
ax2.indicate_inset_zoom(axins2)

ax2.semilogy(A.diagonal())
ax2.set_ylabel('Diagonal')

tight_layout()
savefig('This_Matrix.v4.png')