import numpy as np
import os
import sys
from glm_hmm_utils import *
from neurodatatypes import Session


DPATH = '/u/project/churchland/mmelin/Widefield'
TOL = 10**-4

def crossvalidate(savepath, index):
    #get the params from the npy file
    print(type(index))
    paramspath = os.path.join(savepath,'params.npy')
    params = np.load(paramspath, allow_pickle=True)
    params = params[index,:]
    print(params)

    #[animal, mthd, N_iter, N_state, alpha, sigma]
    animals = params[0]
    sessions = []
    if not isinstance(animals, (list, tuple)):
        animals = [animals]
    for animal in animals:
        sessions_from_one_animal = Session.get_sessions(DPATH, animal, max_nochoice=20, modality=2, min_trials=50, discrim_min=.5, discrim_max=1, assisted_cutoff=.9, singlespout_cutoff=.05)        
        sessions.extend(sessions_from_one_animal)

    glm1 = GlmHmm(sessions,
                  'target_rate',
                  input_terms_list=['coherence','bias'])

    #run the crossvalidation
    log_likelihood = glm1.k_fold_crossval(10, #folds of cross validation
                                          10, #initializations per fold 
                                          n_cpu=10, #ncpus to use for the initializations
                                          scoring='bits',
                                          mthd=params[1],
                                          N_states=int(params[3]),
                                          N_iters=int(params[2]),
                                          tol=TOL,
                                          alpha=int(params[4]),
                                          sigma=float(params[5]),
                                          silent = True)[1]

    #save the output to a file
    savepath = os.path.join(savepath, 'outputs') 
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    savepath = os.path.join(savepath,'{}.npy'.format(index))
    with open(savepath, 'wb') as f:
        np.save(f, np.array(log_likelihood))

if __name__ == '__main__':
    savepath = sys.argv[1] #path for this particular run
    print(sys.argv[2])
    index = int(sys.argv[2]) - 1 # UGE array submission does not allow job ID to start at zero
    crossvalidate(savepath, index)
