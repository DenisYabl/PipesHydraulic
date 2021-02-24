import numpy as np
from scipy.interpolate import interp2d
from itertools import product


def create_lazy_spline_cache_f_wrapper(f, dx=0.1, dy=0.1, half_nx=4, half_ny=4, kind='cubic'):
    nx = 2*half_nx
    ny = 2*half_ny
    grid_step_x = dx*(nx-1)
    grid_step_y = dy*(ny-1)
    covered = set()
    surface_pieces = dict()
    f_on_grid_cache = dict()

    def get_nearest_grid_knot_idxs(x, y):
        x_knot_idx = round(x / grid_step_x)
        y_knot_idx = round(y / grid_step_y)
        return x_knot_idx, y_knot_idx

    def build_spline_surface_over_knot(x_knot_idx, y_knot_idx):
        center_i = x_knot_idx * (nx-1)
        center_j = y_knot_idx * (ny-1)
        ld_corner_i = center_i - half_nx
        ld_corner_j = center_j - half_ny
        ru_corner_i = center_i + half_nx
        ru_corner_j = center_j + half_ny
        i_range = list(range(ld_corner_i, ru_corner_i+1))
        j_range = list(range(ld_corner_j, ru_corner_j+1))
        assert len(i_range) == nx+1
        assert len(j_range) == ny+1
        zs = np.zeros((nx+1, ny+1))
        for (i, j) in product(i_range, j_range):
            if not (i, j) in f_on_grid_cache:
                z = f(i * dx, j * dy)
                f_on_grid_cache[(i, j)] = z
            else:
                z = f_on_grid_cache[(i, j)]
            ii = i - ld_corner_i
            jj = j - ld_corner_j
            zs[ii, jj] = z
        xmin, xmax = i_range[0] * dx, i_range[-1] * dx
        ymin, ymax = j_range[0] * dy, j_range[-1] * dy
        xs = np.linspace(start=xmin, stop=xmax, num=nx+1)
        ys = np.linspace(start=ymin, stop=ymax, num=ny+1)
        zs = np.transpose(zs)
        surface_piece = interp2d(x=xs, y=ys, z=zs, kind=kind)
        return surface_piece

    def f_caller(x, y):
        x_knot_idx, y_knot_idx = get_nearest_grid_knot_idxs(x, y)
        if not (x_knot_idx, y_knot_idx) in covered:
            surface_piece_spline = build_spline_surface_over_knot(x_knot_idx, y_knot_idx)
            covered.add((x_knot_idx, y_knot_idx))
            surface_pieces[(x_knot_idx, y_knot_idx)] = surface_piece_spline
        else:
            surface_piece_spline = surface_pieces[(x_knot_idx, y_knot_idx)]
        z = surface_piece_spline(x, y)
        return z

    return f_caller
