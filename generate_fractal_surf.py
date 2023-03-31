"""
%*************************************************************************%
% All rights reserved (C) to the authors: Mahmoud SHAQFA, Gary CHOI,      %
% Guillaume ANCIAUX & and Katrin BEYER                                    %
%                                                                         %
% M. Shaqfa Contact:                                                      %
% Department of Mechanical Engineering, Massachusetts Institute of        %
% Technology (MIT)                                                        %
% Cambridge, MA, USA                                                      %
%               Email: mshaqfa@mit.edu                                    %
%                                                                         %
% G. Choi Contact:                                                        %
% Department of Mathematics, Massachusetts Institute of Technology (MIT)  %
% Cambridge, MA, USA                                                      %
%               Email: ptchoi@mit.edu                                     %
%                                                                         %
% G. Anciaux Contact:                                                     %
%               Email: guillaume.anciaux@epfl.ch                          %
%                                                                         %
% K. Beyer Contact:                                                       %
%               Email: katrin.beyer@epfl.ch                               %
%*************************************************************************%
% This code includes implementations for:                                 %
%				- Disk harmonics expansion for parametric surfaces        %
% This code is part of the paper: "Disk Harmonics for Analysing Curved    %
% and Flat Self-affine Rough Surfaces and the Topological                 %
% Reconstruction of Open Surfaces"                                        %
%                                                                         %
%*************************************************************************%
% This library is free software; you can redistribute it and/or modify	  %
% it under the terms of the GNU Lesser General Public License as published%
% by the Free Software Foundation; either version 2.1 of the License, or  %
% (at your option) any later version.                                     %
%                                                                         %
% This library is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU Lesser General Public License for more details.        	  %
% You should have received a copy of the GNU Lesser General Public License%
% along with this library; if not, write to the Free Software Foundation, %
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA       %
%*************************************************************************%
% Author of this file: Mahmoud S. Shaqfa, Gary Choi

@author: mahmoudshaqfa

Citing
Tamaas is the result of a science research project. To give proper credit to Tamaas and the researchers who have developed the numerical methods that it implements, please cite Tamaas as:

Frérot , L., Anciaux, G., Rey, V., Pham-Ba, S., & Molinari, J.-F. Tamaas: a library for elastic-plastic contact of periodic rough surfaces. Journal of Open Source Software, 5(51), 2121 (2020). doi:10.21105/joss.02121

If you use the elastic-plastic contact capabilities of Tamaas, please cite:

Frérot, L., Bonnet, M., Molinari, J.-F. & Anciaux, G. A Fourier-accelerated volume integral method for elastoplastic contact. Computer Methods in Applied Mechanics and Engineering 351, 951–976 (2019) doi:10.1016/j.cma.2019.04.006.

If you use the adhesive contact capabilities of Tamaas, please cite:

Rey, V., Anciaux, G. & Molinari, J.-F. Normal adhesive contact on rough surfaces: efficient algorithm for FFT-based BEM resolution. Comput Mech 1–13 (2017) doi:10.1007/s00466-017-1392-5.
"""

import numpy as np
import tamaas as tm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pyvista as pv
from mpl_toolkits.mplot3d import Axes3D


# Helper functions
def make_unit_disk(surface):
    _temp = np.arange(surface.shape[0])
    _xx, _yy = np.meshgrid(_temp, _temp)
    scale = max(surface.shape)
    _xx, _yy, _zz = _xx.ravel(), _yy.ravel(), surface.ravel()

    # Normalize size and location
    _xx = _xx/float(scale)
    _yy = _yy/float(scale)
    if max(_zz) != 0.:
        _zz = 0.07*_zz/max(_zz)
    _xx -= 0.5
    _yy -= 0.5

    # cut a unit disk for mapping
    _temp_radials = np.sqrt(_xx**2 + _yy**2)
    _dels = np.where(_temp_radials > 0.5)
    _xx = np.delete(_xx, _dels)
    _yy = np.delete(_yy, _dels)
    _zz = np.delete(_zz, _dels)
    return _xx, _yy, _zz


def stereographic_hemisphere_projection_from_disk(xd, yd, zd):
    _xx, _yy, _zz = xd, yd, zd

    # Scale ~=0.9999999 (1.0) for a hemisphere
    _scale = np.sqrt(1-np.cos(np.pi/2)/(1+np.cos(np.pi/2)))
    _xx *= _scale
    _yy *= _scale

    # Project on a hemisphere via sterographic projection (angle preserving)
    x = 2*_xx/(1+_xx**2+_yy**2)
    y = 2*_yy/(1+_xx**2+_yy**2)
    z = -(1-_xx**2-_yy**2)/(1+_xx**2+_yy**2)

    # Add the noise along the normal direction
    x = np.multiply(x, (1+_zz))
    y = np.multiply(y, (1+_zz))
    z = np.multiply(z, (1+_zz))
    return x, y, z, _zz


def stereographic_hemisphere_projection(surface):
    _xx, _yy, _zz = make_unit_disk(surface)

    # Scale ~=0.9999999 (1.0) for a hemisphere
    _scale = np.sqrt(1-np.cos(np.pi/2)/(1+np.cos(np.pi/2)))
    _xx *= _scale
    _yy *= _scale

    # Project on a hemisphere via sterographic projection (angle preserving)
    x = 2*_xx/(1+_xx**2+_yy**2)
    y = 2*_yy/(1+_xx**2+_yy**2)
    z = -(1-_xx**2-_yy**2)/(1+_xx**2+_yy**2)

    # Add the noise along the normal direction
    x = np.multiply(x, (1+_zz))
    y = np.multiply(y, (1+_zz))
    z = np.multiply(z, (1+_zz))
    return x, y, z, _zz


save_hemisphere = False
# Create an Isotropic power law
spectrum = tm.Isopowerlaw2D()

# Set spectrum parameters
n = 2**8
spectrum.q0 = 0
spectrum.q1 = 4
spectrum.q2 = 2**8
spectrum.hurst = 0.95
RMS = 10.
seed_no = 9700

RMS_computed = spectrum.rmsHeights()
RMS_slope = spectrum.rmsSlopes()
print("RMS of the surface is: {}, and the RMS slope is: {}".format(RMS_computed, RMS_slope))

# Generate the surface
# generator = tm.SurfaceGeneratorFilter2D()     # Generate surface with randomized phase and amplitude
generator = tm.SurfaceGeneratorRandomPhase2D()  # Generate surface with only random phase
if seed_no is not None:
    generator.random_seed = seed_no
generator.setSizes([n, n])
generator.setFilter(spectrum)
surface = generator.buildSurface()

u = np.linspace(0, n-1, endpoint=True, num=n)
v = np.linspace(0, n-1, endpoint=True, num=n)
u, v = np.meshgrid(u, v)

# Update RMS of the surface
arithmetic_RMS = np.sqrt(np.sum(surface - np.mean(surface))**2/n**2)
print("The arithmetic RMS is (Before modification): {}".format(arithmetic_RMS))
rms_multiplier = RMS/RMS_computed
print(RMS_computed)
surface *= rms_multiplier
arithmetic_RMS = np.sqrt(np.sum(surface - np.mean(surface))**2/n**2)
print("The arithmetic RMS is (After modification): {}".format(arithmetic_RMS))

# Compute RMS statistics from Tamaas
RMS_Tamaas = tm.Statistics2D.computeRMSHeights(surface)
print("RMS heights computed from Tamaas: {}".format(RMS_Tamaas))

# Plot surface
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(u, v, surface)
plt.show()

# Plot point cloud
points = np.c_[u.reshape(-1), v.reshape(-1), surface.reshape(-1)]
cloud = pv.PolyData(points)
surf = cloud.delaunay_2d()
#surf.plot(cpos="xy", show_edges=True)

# cloud.plot(point_size=15)
pv.save_meshio("rec_artificial_H_" + str(spectrum.hurst) + ".stl", surf)
print("Saved the surface")

# Save a unit disk
xd, yd, zd = make_unit_disk(surface)
points = np.c_[xd, yd, zd]
cloud = pv.PolyData(points)
surf = cloud.delaunay_2d()
surf.plot(cpos="xy", show_edges=True)
pv.save_meshio("disk_rec_artificial_H_" + str(spectrum.hurst) + ".stl", surf)
print("Saved the unit disk surface")

# Save a hemisphere
if save_hemisphere:
    xh, yh, zh, _ = stereographic_hemisphere_projection_from_disk(
        points[:, 0], points[:, 1], points[:, 2])
    points = np.c_[xh, yh, zh]
    cloud = pv.PolyData(points)
    #cloud = pv.PolyData(points)
    surf = cloud.delaunay_2d()
    surf.plot(cpos="xy", show_edges=True)
    #pv.save_meshio("hemisphere_rec_artificial_H_" + str(spectrum.hurst) + ".stl", surf)
    print("Saved the unit hemispherical surface")
