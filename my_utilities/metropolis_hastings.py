import numpy as np

def _apply_pbc(x, X_MIN, X_MAX):
    return ((x - X_MIN) % (X_MAX - X_MIN)) + X_MIN

def gen_random_sample(density,
                      N_SAMPLES,
                      N_DIM,
                      START_POINT,
                      RANDOM_STEP_SIGMA,
                      X_MIN = np.array([-np.pi, -np.pi]),
                      X_MAX = np.array([np.pi, np.pi])
                      ):
    
    p = np.array(START_POINT)
    prob_curr = density(p)

    sample = np.empty((N_SAMPLES, N_DIM))

    for i in range(N_SAMPLES):
        # find new random point using
        # normal distribution with center in last point
        p_proposition = np.random.normal(p, RANDOM_STEP_SIGMA)
        # apply periodic boundary condition
        p_proposition = _apply_pbc(p_proposition, X_MIN, X_MAX)
        # what is the probability density for the proposed point?
        prob_new = density(p_proposition)
        
        acceptance_ratio = prob_new / prob_curr
        
#         print(f"{p = }\n{p_proposition = }\n{prob_curr}\n{prob_new}\n{acceptance_ratio}\n")
#         input()
#         print()
        
        if acceptance_ratio > np.random.rand():
            # accept
            p = p_proposition
            prob_curr = prob_new

        sample[i] = p
        
    return sample