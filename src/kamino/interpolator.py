from scipy.interpolate import RegularGridInterpolator
from collections.abc import Callable
import numpy as np
import itertools
from tqdm import tqdm
import os

def make_interpolator(function: Callable, grid: list[np.typing.NDArray[np.float64]], logspace: list[bool] | None=None, interpolator_data_file: str='') -> Callable:

    if not os.path.exists(interpolator_data_file):

        table_shape = [len(arr) for arr in grid]
        indexes = [range(l) for l in table_shape]

        if logspace is None:
            logspace = [False] * len(grid)
        elif len(logspace) != len(grid):
            raise ValueError(f"logspace length ({len(logspace)}) must match grid length ({len(grid)})")
        
        first_args = [arr[0] for arr in grid]
        sample_output = np.asarray(function(*first_args))

        full_shape = table_shape + list(sample_output.shape)
        interpolator_values = np.zeros(full_shape)

        for index, args in tqdm(zip(itertools.product(*indexes), itertools.product(*grid)), total=int(np.prod(table_shape))):
            interpolator_values[index] = function(*args)

        new_grid = [(np.log10(arr) if log else arr) for log, arr in zip(logspace, grid)]

        if interpolator_data_file != '':
            np.savez_compressed(
                interpolator_data_file,
                grid_data=np.array(new_grid, dtype=object), 
                values=interpolator_values,
                logspace=np.array(logspace)
            )

    else:

        data = np.load(interpolator_data_file, allow_pickle=True)
        new_grid = data['grid_data'].tolist()
        interpolator_values = data['values']
        logspace = data['logspace'].tolist()

    RGI = RegularGridInterpolator(new_grid, interpolator_values)

    def interpolator(*args):
        args = [(np.log10(arg) if log else arg) for arg, log in zip(args, logspace)]
        return RGI(args)

    return interpolator

if __name__ == '__main__':

    def dummy_function(a1=0, a2=0, a3=0):
        return a1 * a2, a1 + a2

    f = make_interpolator(dummy_function, [np.linspace(0.1, 1, num=10), np.linspace(0, 2, num=5)], [True, False])

    print(f(0.5, 1.6))
