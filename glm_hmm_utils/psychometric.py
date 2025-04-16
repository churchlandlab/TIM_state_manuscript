import numpy as np
from scipy.optimize import minimize

def compute_p_right(coherence, choices):
    from statsmodels.stats.proportion import proportion_confint

    coherence_set = np.unique(coherence)
    p_right = np.zeros_like(coherence_set)
    ci_right = np.zeros((len(coherence_set),2))

    for i,coh in enumerate(coherence_set):
        coh_choices = np.squeeze(choices[coherence == coh])
        n_trials = coh_choices.size
        n_right = np.sum(coh_choices == 1)
        p_right[i] = n_right / n_trials
        ci_right[i] = proportion_confint(n_right, n_trials, method='wilson')

    correct_side = np.sign(coherence)
    correct_side[correct_side == -1] = 0
    total_p_correct = np.sum(correct_side == choices) / np.size(correct_side)
    print(f'Performance was: {total_p_correct}')
    return p_right,ci_right

def cumulative_gaussian(alpha,beta,gamma,lmbda, X):
    '''
    Evaluate the cumulative gaussian psychometric function.
       alpha is the bias (left or right)
       beta is the stepness
       gamma is the left handside offset
       lmbda is the right handside offset
      
    Adapted from the Palamedes toolbox 
    Joao Couto - Jan 2022    
    '''
        
    from scipy.special import erfc # import the complementary error function
    return  gamma + (1 - gamma - lmbda)*0.5*erfc(-beta*(X-alpha)/np.sqrt(2))+1e-9    


def neg_log_likelihood_error(func, parameters, X, Y):
    '''
    Compute the log likelihood

    'func' is the (psychometric) function 
    'parameters' are the input parameters to 'func'
    'Y' is the binary response (correct = 1; incorrect=0)
    '''

    pX = func(*parameters, X)*0.99 + 0.005  # the predicted performance for X from the PMF
    # epsilon to prevent error in log(0)
    val = np.nansum(Y*np.log(pX) + (1-Y)*np.log(1-pX))
    return -1*val
