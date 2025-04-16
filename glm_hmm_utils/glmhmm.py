# -*- coding: utf-8 -*-
#from .datahandling import * #TODO: delete this file eventually? it is becoming deprecated
from tkinter import W
import numpy as np
import matplotlib.pyplot as plt
import ssm
import os
import pickle
import glob
from datetime import datetime
import multiprocess as mp

class GlmHmm():
    """
    class GlmHmm is used to create and train GLM-HMM's on the Musall/Sun/Gluf dataset. It 
    is also built on top of the SSM class.
    
    By Max Melin, mmelin@g.ucla.edu
    
    TODO: docstrings for all functions
    TODO: Configure setup.py for new modules
    """
    def __init__(self, Sessions, init_mthd, input_terms_list, obs_dim = 1,
                 num_categories = 2, input_dim = 2):
        """
        Initialize the GlmHmm class.
        Sessions is a list of Session objects.
        """
        self.inpts = None
        self.outputs = None
        self.obs_dim = obs_dim
        self.num_categories = num_categories # should always be 2 for binary left/right choice
        self.input_terms_list = input_terms_list
        self.input_dim = len(input_terms_list)
        self.nummodels = 0
        self.models = []
        self.train_log_likelihoods = []
        self.training_params = []
        self.init_mthd = init_mthd
        if Sessions is not None:
            self._generate_model_inputs(Sessions) #generates self.inpts and self.outputs
            self.session_dates = [sess.date for sess in Sessions]
            self.session_animals = [sess.animal for sess in Sessions]

    def create_untrained_model(self, mthd, N_states, alpha=2, sigma=2):
            
        if mthd == 'mle':
            ssmobject = ssm.HMM(N_states, self.obs_dim, self.input_dim, observations="input_driven_obs", 
                           observation_kwargs=dict(C=self.num_categories), transitions="standard")

        elif mthd == 'map':
            ssmobject = ssm.HMM(N_states, self.obs_dim, self.input_dim, observations="input_driven_obs", 
                         observation_kwargs=dict(C=self.num_categories,prior_sigma=sigma),
                         transitions="sticky", transition_kwargs=dict(alpha=alpha,kappa=0))
            
        
        self.models.append(ssmobject)
        self.nummodels += 1
            
    def return_best_model(self, scoring, inpts, outputs, train_outputs=None, silent=True):
        all_init_lls = []
        for mdl in self.models:
            if scoring == 'll':
                all_init_lls.append(mdl.log_likelihood(outputs, inputs=inpts))
            elif scoring == 'bits':
                assert train_outputs is not None, 'Must specify train_outputs for bits scoring'
                session_lengths = [arr.size for arr in outputs]
                n_trials = np.sum(session_lengths)
                ll = mdl.log_likelihood(outputs, inputs=inpts)
                ll_0 = self._calculate_baseline_test_ll(train_outputs, outputs, self.num_categories)
                bits_per_trial = self._calculate_cv_bit_trial(ll, ll_0, n_trials)
                all_init_lls.append(bits_per_trial)
            else:
                raise ValueError('Please specify a valid scoring method')
        all_init_lls = np.array(all_init_lls)
        best_model_index = np.argmax(all_init_lls)

        if not silent:
            print(f'Log likelihoods for initializations are {all_init_lls}')
        print(f'Model index {best_model_index} has the highest log likelihood: {all_init_lls[best_model_index]}')
        return self.models[best_model_index], best_model_index, all_init_lls[best_model_index]
    
    def _clear_saved_models(self):
        self.models = []
        self.train_log_likelihoods = []
        self.training_params = []
        self.nummodels = 0
            
    def train_n_models(self, n_models, n_cpu='all', **train_new_model_kwargs):
        #train n models in parallel. important for random initialization
        if n_cpu == 'all':
            n_cpu = mp.cpu_count()

        def mpfunc(iterator):
            np.random.seed() #linux forks processes, so need to reseed np.random or the results from each worker will be the same
            return self._train_new_model(**train_new_model_kwargs)
        
        if n_cpu == 1:
            res = list(map(mpfunc, range(n_models)))
        else:
            with mp.Pool(n_cpu) as p:
                res = p.map(mpfunc, range(n_models)) # pass by reference
        trained_models, train_log_likes = zip(*res)

        self.models.extend(trained_models) #save the trained ssm object
        self.train_log_likelihoods.extend(train_log_likes) #save the trained ssm object
        return trained_models

    def train_new_model(self, mthd, **kwargs):
        if self.nummodels == 0:
            print('\nTraining the first model for this instance.')
        else:
            print('\ntraining model number {} for this instance.'.format(self.nummodels+1))

        trained_model, ll = self._train_new_model(mthd, **kwargs)
        self.models.append(trained_model)
        self.log_likelihoods.append(ll)

        self.nummodels += 1
        return trained_model

    def _train_new_model(self, mthd=None, inpts = None, outputs = None, N_iters = 200, N_states = 3, tol = 10**-6 ,  alpha = 2, sigma = 2, silent = False):
       
        if inpts or outputs is not None:
            print('\nUsing a subselection of all inputs')
            #numsess = len(inpts)
        else:
            #print('\nUsing all specified dates from when the glmhmm object was created: ' + str(self.dates) + '\n')
            print('\nUsing all specified dates from when the glmhmm object was created. \n')
            #numsess = len(self.dates) 
            inpts = self.inpts
            outputs = self.outputs
            
        assert mthd is not None, 'Please specify a training method'
        
        if mthd == 'mle':
            ssmobject = ssm.HMM(N_states, self.obs_dim, self.input_dim, observations="input_driven_obs", 
                           observation_kwargs=dict(C=self.num_categories), transitions="standard")

            log_likelihood = ssmobject.fit(outputs, inputs=inpts, method='em', num_iters=N_iters, tolerance=tol) 

        elif mthd == 'map':
            ssmobject = ssm.HMM(N_states, self.obs_dim, self.input_dim, observations="input_driven_obs", 
                         observation_kwargs=dict(C=self.num_categories,prior_sigma=sigma),
                         transitions="sticky", transition_kwargs=dict(alpha=alpha,kappa=0))
            
            log_likelihood = ssmobject.fit(outputs, inputs=inpts, method='em', num_iters=N_iters, tolerance=tol)
        else:
            raise ValueError('Please specify a valid training method')
        
        training_data_dict = {'alpha':alpha,
                              'sigma':sigma,
                              'mthd':mthd,
                              'N_iters':N_iters,
                              'N_states':N_states,
                              'tol':tol}
        
        #self.training_params.append(training_data_dict) # FIXME: move to outer function
    
        if not silent: # Plotting
            fig = plt.figure(figsize=(4, 3), dpi=80, facecolor='w', edgecolor='k')
            plt.plot(log_likelihood)
            plt.title('Log Likelihood ({} method)'.format(mthd))
            plt.xlabel("EM Iteration")
            plt.ylabel("Log Probability")
            plt.show()
        
        return ssmobject, log_likelihood[-1] #returns the ssm object, and ll at the end of training

    def train_initialized_model(self, mthd, inpts = None, outputs = None, N_iters = 200, N_states = 3, tol = 10**-9 ,  alpha = 2, sigma = 2, silent = False):
       
        if self.nummodels == 0:
            print('\nThere are no models to initialize from, exiting now.')
            return
        else:
            print('\nInitializing parameters from previous model and training model number {} for this instance.'.format(self.nummodels+1))  
        
        if inpts or outputs is not None:
            print('\nUsing a subselection of all inputs')
            numsess = len(inpts)
        else:
            print('\nUsing all specified dates from when the glmhmm object was created: ' + str(self.dates) + '\n')
            numsess = len(self.dates) 
            inpts = self.inpts
            outputs = self.outputs
            
        if mthd == 'mle':
            ssmobject = ssm.HMM(N_states, self.obs_dim, self.input_dim, observations="input_driven_obs", 
                           observation_kwargs=dict(C=self.num_categories), transitions="standard")


        elif mthd == 'map':
            ssmobject = ssm.HMM(N_states, self.obs_dim, self.input_dim, observations="input_driven_obs", 
                         observation_kwargs=dict(C=self.num_categories,prior_sigma=sigma),
                         transitions="sticky", transition_kwargs=dict(alpha=alpha,kappa=0))
            
        ssmobject.params = self.models[-1].params
        log_likelihood = ssmobject.fit(outputs, inputs=inpts, method='em', num_iters=N_iters, tolerance=tol, initialize=False) # must set initialization to false to use the parameters we give
        
        self.models.append(ssmobject) #save the trained ssm object
        
        training_data_dict = {'alpha':alpha,
                              'sigma':sigma,
                              'mthd':mthd,
                              'N_iters':N_iters,
                              'N_states':N_states,
                              'tol':tol}
        
        self.training_params.append(training_data_dict) #save training parameters for that object
        self.nummodels += 1
    
        if not silent: # Plotting
            fig = plt.figure(figsize=(4, 3), dpi=80, facecolor='w', edgecolor='k')
            plt.plot(log_likelihood)
            plt.title('Log Likelihood ({} method)'.format(mthd))
            plt.xlabel("EM Iteration")
            plt.ylabel("Log Probability")
            plt.show()
        
        return ssmobject

    def k_fold_crossval(self, n_folds, n_initializations_per_fold, scoring='bits', **training_kwargs):
        from sklearn.model_selection import KFold
        if n_folds is None:
            n_folds = len(self.dates) #hold out one session for each fold of cross-validation
        session_inds = np.arange(len(self.inpts))
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=None)
        
        inpts = np.array(self.inpts, dtype=object) # need to cast as ndarray for indexing, will convert back to list later
        outputs = np.array(self.outputs, dtype=object)
        
        log_likelihoods = np.array([])
        
        self._clear_saved_models()
        for train_index, test_index in kf.split(session_inds):
            print("TRAIN:", train_index, "TEST:", test_index)
            train_inpt = inpts[train_index].tolist()
            train_output = outputs[train_index].tolist()
            test_inpt = inpts[test_index].tolist()
            test_output = outputs[test_index].tolist()
            
            #mdl = self.train_new_model(mthd, inpts = train_inpt, outputs = train_output, 
            #                       N_iters = N_iters, N_states = N_states, tol = tol, alpha = alpha,
            #                       sigma = sigma, silent = silent)
            
            self.train_n_models(n_initializations_per_fold,
                                inpts=train_inpt,
                                outputs=train_output,
                                **training_kwargs)
            mdl, best_mdl_index, ll = self.return_best_model(scoring,
                                                             inpts=test_inpt,
                                                             outputs=test_output,
                                                             train_outputs=train_output)
            log_likelihoods = np.append(log_likelihoods, ll)
            self._clear_saved_models() #wipe the models that were just trained before moving to the next fold
            
        return np.mean(log_likelihoods), log_likelihoods

    def _calculate_cv_bit_trial(self, ll_model, ll_0, n_trials):
        '''From Ashwood et al.'''
        cv_bit_trial = ((ll_model - ll_0) / n_trials) / np.log(2)
        return cv_bit_trial
    
    def _calculate_baseline_test_ll(self, train_y, test_y, C):
        """
        from Ashwood et al.
        Calculate baseline loglikelihood for CV bit/trial calculation.  This is
        log(p(y|p0)) = n_right(log(p0)) + (n_total-n_right)log(1-p0), where p0
        is the proportion of trials
        in which the animal went right in the training set and n_right is the
        number of trials in which the animal went right in the test set
        :param train_y
        :param test_y
        :return: baseline loglikelihood for CV bit/trial calculation
        """
        train_y = np.vstack(train_y)
        test_y = np.vstack(test_y) # need to stack the sessions for the calculation

        _, train_class_totals = np.unique(train_y, return_counts=True)
        train_class_probs = train_class_totals / train_y.shape[0]
        _, test_class_totals = np.unique(test_y, return_counts=True)
        ll0 = 0
        for c in range(C):
            ll0 += test_class_totals[c] * np.log(train_class_probs[c])
        return ll0


    def return_ordered_weights(self, modelindex = None):
        if modelindex == None:
            modelindex = len(self.models)-1
        unordered_weights = self.models[modelindex].observations.params.squeeze()
        if not hasattr(self, 'label_inds'):
            print('\nUser must first generate labels for the states')
            self.user_label_states(modelindex=modelindex, silent=False)
        ordered_weights = unordered_weights[self.label_inds,:]
        print('\nThese are the ordered weights, make sure they look right:')
        print(ordered_weights)
        return ordered_weights
    
    def return_ordered_states(self, modelindex = None):
        if modelindex == None:
            modelindex = len(self.models) - 1
        posterior_probs_unordered = [self.models[modelindex].expected_states(data=data, input=inpt)[0] for data, inpt in zip(self.outputs, self.inpts)] 
        if not hasattr(self, 'label_inds'):
            print('\nUser must first generate labels for the states')
            self.user_label_states(modelindex = modelindex)

        posterior_probs_ordered = []
        for session in posterior_probs_unordered:
            session = session[:,self.label_inds]
            posterior_probs_ordered.append(session)
        return posterior_probs_ordered
        
    def user_label_states(self, silent=True, modelindex = None):
        if modelindex == None:
            modelindex = len(self.models)-1 #plots the last model trained by default
        unordered_weights = self.models[modelindex].observations.params.squeeze().T
        if not silent:
            plt.plot(unordered_weights,'-o')
            plt.show()
        entered_list = input('Please input the indices of engaged, left biased, and right biased, starting at 0, separated by a space.').split()
        #TODO: assert that entered_list is the same length as num_states for the model
        entered_list = list(map(int,entered_list))
        #print(f'\nEngaged index is {entered_list[0]}, left bias index is {entered_list[1]}, right bias index is {entered_list[2]}') 
        
        self.label_inds = entered_list  #[engaged_ind, lbias_ind, rbias_ind]

    def _generate_model_inputs(self,Sessions):
        # Get and format data from multiple sessions
        allchoices = [] # a list of numpy arrays, where the list is N sessions long
        inpts = []# a list of numpy arrays, where the list is N sessions long
        masks = [] 
        for i,sess in enumerate(Sessions):
            inpt, choice, mask = self.__format_session(sess)
            allchoices.append(choice)
            inpts.append(inpt)
            masks.append(mask)
        
        self.inpts = inpts
        self.outputs = allchoices
        self.masks = masks
        self.dates = [sess.date for sess in Sessions]

    def __format_session(self,Session): #TODO: move to datahandling.py to replace deprecated function once this is known to work. 
        
        

        if self.init_mthd == 'target_rate':
            #this is using the target rate for left and right stimuli
            stimside = Session.data.stimside.values
            targstim = Session.data.targstim.values
            diststim = Session.data.diststim.values

            lstim = np.zeros(len(stimside))
            rstim = np.zeros(len(stimside))
            lstim[stimside==0] = targstim[stimside==0]
            lstim[stimside==1] = diststim[stimside==1]
            rstim[stimside==0] = diststim[stimside==0]
            rstim[stimside==1] = targstim[stimside==1]

            coherence = rstim - lstim

        elif self.init_mthd == 'actual_rate':
            # this is using the actual rate, rather than the target rate 
            allleftclicks = np.empty(0)
            allrightclicks = np.empty(0)
            allleftflashes = np.empty(0)
            allrightflashes = np.empty(0)
        
            for i in range(num_trials): 
                numleftclicks = np.size(Session.stimtimes[i][0]) #TODO: this may need to be changed to what's reflected in NeuroDataTypes?
                numrightclicks = np.size(Session.stimtimes[i][1])
                #numleftflashes = np.size(Session.data.stimstruct[i][2]) not using visual data yet
                #numrightflashes = np.size(Session.data.stimstruct[i][3])

                allleftclicks = np.append(allleftclicks,numleftclicks)
                allrightclicks = np.append(allrightclicks,numrightclicks)
                #allleftflashes = np.append(allleftflashes,numleftflashes)
                #allrightflashes = np.append(allrightflashes, numrightflashes)
            
            coherence = allrightclicks - allleftclicks
            coherence = coherence[trial_mask]

        else:
            ValueError('Please input a valid init_method.')

        coherence = -coherence #sign flip so that weights come out positive


        choice = Session.data.choice.values
        discrimination = Session.data.discrimination.values
        num_trials = len(choice)

        
        inpt = np.empty([num_trials,self.input_dim])
        for i, input_term in enumerate(self.input_terms_list):
            if input_term == 'coherence':
                inpt[:,i] = coherence
                coherence_index = i
            elif input_term == 'bias':
                inpt[:,i] = np.negative(np.ones([num_trials]))
            elif input_term == 'previous_choice':
                inpt[:,i] = np.concatenate([[np.nan], choice[:-1]])
            elif input_term == 'wsls':
                inpt[:,i] = np.concatenate([[np.nan], ~stimside[:-1].astype(np.bool)]) # previous stimside = win stay, lose switch
        
        choice_mask = ~np.isnan(choice) # remove trials with no choice and with opto
        disc_mask = discrimination == 1
        num_nochoice = np.sum(np.isnan(choice))
        invalid_inputs = np.any(np.isnan(inpt), axis=1)

        print(f'\n{num_trials} trials total. Removing {num_nochoice} trials with no choice.')
        print(f'Removing {int(num_trials - np.sum(discrimination))} trials that are set to detection.')
        #trial_mask = np.logical_and(choice_mask, disc_mask)
        trial_mask = np.logical_and.reduce((choice_mask, disc_mask, ~invalid_inputs))
        choice = choice[trial_mask]
        choice = np.int32(choice) # ssm module doesn't like floats
        choice = np.expand_dims(choice, axis=1)
        inpt = inpt[trial_mask,:]
        num_trials = choice.size  #number of trials per session with choice

        inpt[:,coherence_index] = inpt[:,coherence_index]/ np.max(inpt[:,coherence_index]) #normalize after discarding bad trials
        #from sklearn.preprocessing import scale
        #coherence = scale(coherence)
        return inpt, choice, trial_mask
        
    def plot_transition_matrix(self):
        raise NotImplementedError()
        #TODO: make this code work, it's just been directly copied over so far
        fig = plt.figure(figsize=(5, 2.5), dpi=80, facecolor='w', edgecolor='k')
        recovered_trans_mat = np.exp(glmhmm.transitions.log_Ps)
        plt.imshow(recovered_trans_mat, vmin=-0.8, vmax=1, cmap='bone')
        for i in range(recovered_trans_mat.shape[0]):
            for j in range(recovered_trans_mat.shape[1]):
                text = plt.text(j, i, str(np.around(recovered_trans_mat[i, j], decimals=2)), ha="center", va="center",
                                color="k", fontsize=12)
        plt.xlim(-0.5, num_states - 0.5)
        plt.xticks(range(0, num_states), ('1', '2', '3'), fontsize=10)
        plt.yticks(range(0, num_states), ('1', '2', '3'), fontsize=10)
        plt.ylim(num_states - 0.5, -0.5)
        plt.title(str(mthd), fontsize = 15)
        plt.subplots_adjust(0, 0, 1, 1)

    def publication_plot(self):
        raise NotImplementedError()
        #this line will be useful for pdf saving
        fig.savefig(r'C:\Data\churchland\myimage.pdf', format='pdf', dpi=1000, bbox_inches='tight')
        
    def save(self, savedir, modelindex=None):
        import scipy.io as sio
        if modelindex is None:
            modelindex = len(self.models)-1
        #TODO: check the code below and modify as needed. maybe allow to save multiple models to seperate files? Or save the whole large object?
        savename = input('Enter a name for saving data: ')
        if not hasattr(self, 'label_inds'):
            print('\nUser must first generate labels for the states')
            self.user_label_states(modelindex = modelindex)
            print(self.label_inds)

        posterior_probs = [self.models[modelindex].expected_states(data=data, input=inpt)[0] for data, inpt in zip(self.outputs, self.inpts)] 

        #TODO: add saving animal to the dictionary
        savedic = {'glmhmm_params': self.models[modelindex],
                   'masks': self.masks,
                   'state_label_indices': list(np.array(self.label_inds) + 1), # add 1 for matlab indexing
                   'mouse': self.session_animals,
                   'model_training_sessions': self.session_dates,
                   'posterior_probs': posterior_probs}

        sio.savemat(os.path.join(savedir, savename + '.mat'),savedic)
        with open(os.path.join(savedir, savename + '.pickle'), 'wb') as handle:
            #pickle.dump(self.models[modelindex], handle, protocol=pickle.HIGHEST_PROTOCOL)
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)


    def save_batch_run(self, savedir, savename):
        with open(os.path.join(savedir, savename + '.pickle'), 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        
if __name__ == '__main__':
    from neurodatatypes import Session
    DPATH = r'/u/project/churchland/mmelin/Widefield'
    mouse = 'mSM63'
    sessionlist = Session.get_sessions(DPATH, mouse, max_nochoice = 5, modality = 2, min_trials = 200, discrim_min = .77, discrim_max = 1, assisted_cutoff = .95, singlespout_cutoff = .02)
    model = GlmHmm(sessionlist, 'target_rate')
    
        
        