import numpy as np

def log_likelihood(density_function, data):
    return np.sum(np.array([np.log(density_function(point)) for point in data]))