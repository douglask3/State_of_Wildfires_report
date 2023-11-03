import numpy as np

def round_to_sf(numbers, ndigs = 2):
    def sf(number):
        if number == 0:
            return 0  # Handle the special case of 0
        order = int(np.floor(np.log10(np.abs(number))))
        factor = 10 ** (ndigs - 1 - order)
        return round(number * factor) / factor
    return np.array([sf(number) for number in numbers])


