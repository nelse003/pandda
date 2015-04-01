import os, sys, glob, time, re

import numpy

from cctbx import maptbx
from cctbx import crystal

from scitbx.array_family import flex

from Giant.Grid.Utils import get_grid_points_within_distance_cutoff_of_origin, combine_grid_point_and_grid_vectors

class atomic_mask(object):
    def __init__(self, cart_sites, grid_size, unit_cell, max_dist, min_dist):
        """Take a grid and calculate all grid points with a certain distance cutoff of any point in cart_sites"""

        if min_dist: assert max_dist > min_dist, 'Minimum Mask Distance must be smaller than Maximum Mask Distance'
        print('===================================>>>')

        # TODO Check that all of the points lie withing the fake unit cell
        # i.e. cart_sites + max_dist < unit_cell.parameters()

        # Store distances from masking atoms
        self._max_dist = max_dist
        self._min_dist = min_dist

        # Store grid size
        self._grid_size = grid_size
        self._grid_idxr = flex.grid(grid_size)

        # Unit cell
        self._fake_unit_cell = unit_cell

        t1 = time.time()

        # Calculate the masked indices defined by max distance from protein atoms
        self._outer_mask_indices = maptbx.grid_indices_around_sites(unit_cell=unit_cell,
                                    fft_n_real=grid_size, fft_m_real=grid_size,
                                    sites_cart=cart_sites, site_radii=flex.double(cart_sites.size(), max_dist))
        # Calculate a binary list of the selected indices
        outer_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=int)
        [outer_mask_binary.put(idx, 1) for idx in self.outer_mask_indices()]
        self._outer_mask_binary = outer_mask_binary.tolist()
        # Select the masked grid points
        self._outer_mask_points = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._outer_mask_binary[idx]==1]

        t2 = time.time()
        print('> OUTER MASK > Time Taken: {!s} seconds'.format(int(t2-t1)))

        # Calculate the masked indices defined by min distance from protein atoms
        self._inner_mask_indices = maptbx.grid_indices_around_sites(unit_cell=unit_cell,
                                    fft_n_real=grid_size, fft_m_real=grid_size,
                                    sites_cart=cart_sites, site_radii=flex.double(cart_sites.size(), min_dist))
        # Calculate a binary list of the selected indices
        inner_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), dtype=int)
        [inner_mask_binary.put(idx, 1) for idx in self.inner_mask_indices()]
        self._inner_mask_binary = inner_mask_binary.tolist()
        # Select the masked grid points
        self._inner_mask_points = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._inner_mask_binary[idx]==1]

        t3 = time.time()
        print('> INNER MASK > Time Taken: {!s} seconds'.format(int(t3-t2)))

        # Calculate the combination of these masks
        total_mask_binary = numpy.zeros(self._grid_idxr.size_1d(), int)
        [total_mask_binary.put(self._grid_idxr(gp), 1) for gp in self.outer_mask()]
        [total_mask_binary.put(self._grid_idxr(gp), 0) for gp in self.inner_mask()]
        self._total_mask_binary = total_mask_binary.tolist()
        # Select the masked grid indices
        self._total_mask_indices = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._total_mask_binary[idx]==1]
        # Select the masked grid points
        self._total_mask_points = [gp for idx, gp in enumerate(flex.nested_loop(self._grid_size)) if self._total_mask_binary[idx]==1]

        t4 = time.time()
        print('> TOTAL MASK > Time Taken: {!s} seconds'.format(int(t4-t3)))

    def total_mask(self):
        """Return the grid points allowed by the mask - combination of max_dist (allowed) and min_dist (rejected)"""
        return self._total_mask_points
    def outer_mask(self):
        """Get grid points allowed subject to max_dist"""
        return self._outer_mask_points
    def inner_mask(self):
        """Get grid points rejected subject to min_dist"""
        return self._inner_mask_points

    def total_mask_binary(self):
        return self._total_mask_binary
    def outer_mask_binary(self):
        return self._outer_mask_binary
    def inner_mask_binary(self):
        return self._inner_mask_binary

    def total_mask_indices(self):
        return self._total_mask_indices
    def outer_mask_indices(self):
        return self._outer_mask_indices
    def inner_mask_indices(self):
        return self._inner_mask_indices

    def total_size(self):
        """Returns the number of grid points in the mask"""
        return len(self.total_mask())
    def outer_size(self):
        """Returns the number of grid points inside max_dist"""
        return len(self.outer_mask())
    def inner_size(self):
        """Returns the number of grid points inside max_dist"""
        return len(self.inner_mask())

    def extent(self):
        """Returns the minimum and maximum grid points in the mask"""
        return min(self.total_mask()), max(self.total_mask())

    def summary(self):
        return '\n'.join([  '===================================>>>',
                            'Atomic Mask Summary:',
                            'Total Mask Size (1D): {!s}'.format(self.total_size()),
                            'Outer Mask Size (1D): {!s}'.format(self.outer_size()),
                            'Inner Mask Size (1D): {!s}'.format(self.inner_size()),
                            'Masked Grid Min/Max: {!s}'.format(self.extent())
                        ])

class spherical_mask(object):
    def __init__(self, grid_spacing, distance_cutoff, grid_jump=None):
        """Sphere used to mask grid points within a certain distance of a point"""

        self._mask = get_grid_points_within_distance_cutoff_of_origin(grid_spacing=grid_spacing, distance_cutoff=distance_cutoff)
        self._radius = distance_cutoff
        self._buffer = max(max(self._mask))
        self._grid_spacing = grid_spacing
        if grid_jump:
            self._grid_jump = int(grid_jump)
        else:
            self._grid_jump = iceil(self._buffer)
        if self._grid_jump == 0:
            self._grid_jump = 1

    def mask(self):
        return self._mask
    def size(self):
        return len(self.mask())
    def buffer_size(self):
        return self._buffer
    def grid_spacing(self):
        return self._grid_spacing
    def grid_jump(self):
        return self._grid_jump
    def radius(self):
        return self._radius
    def volume(self):
        return (4/3.0)*numpy.pi*(self.radius()**3)

    def apply_mask(self, grid_point):
        """Combine a grid point with all of the masking vectors"""
        return combine_grid_point_and_grid_vectors(start_point=grid_point, grid_vectors=self.mask())

    def summary(self):
        return '\n'.join(['===================================>>>',
                          'Local Mask Summary:',
                          'Number of Mask Points:  {!s}'.format(len(self.mask())),
                          'Mask Radius (Cart):     {!s}'.format(self.radius()),
                          'Mask Volume (Cart):     {!s}'.format(round(self.volume(),3)),
                          'Largest Mask Vector:    {!s}'.format(max(self.mask())),
                          'Req. Edge Buffer Zone:  {!s}'.format(self.buffer_size()),
                          'Sampling Fraction (1D): 1/{!s}'.format(self.grid_jump()),
                          'Sampling Fraction (3D): 1/{!s}'.format(self.grid_jump()**3)
                        ])


