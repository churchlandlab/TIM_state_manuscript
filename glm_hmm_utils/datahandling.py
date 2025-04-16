# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 15:21:11 2021
@author: mmelin, mmelin@g.ucla.edu

Dataset specific functions used to format data for the GlmHmm object

Run this module as __main__ to get all session info for a mouse and place into a text file.
"""
import numpy as np
import scipy.io as sio
import os
import glob
from datetime import datetime

#TODO: update docstrings
    
def createSessionTxtFile(dpath,mouse):
    """
    Generate a .txt file with information about all sessions for a given subject

    Parameters
    ----------
    dpath : TYPE
        Path to the data.
    mouse : TYPE
        The desired animal ID found in dpath.

    Returns
    -------
    session_dates : TYPE
        list of session dates for that animal.
    modalities : TYPE
        list of modalities for each session.

    """
    datadir = os.path.join(dpath,mouse,'SpatialDisc')
    savefile = os.path.join(dpath,mouse,'SesionData.txt')
    session_dates = [name for name in os.listdir(datadir) if os.path.isdir(os.path.join(datadir, name))] #get directories that contain behavioral session data
    session_dates = sorted(session_dates,key=lambda date: datetime.strptime(date[:11], "%d-%b-%Y")) #parse and sort dates
    numsess = len(session_dates)
    file = open(savefile, 'w')
    file.writelines(['Session info for ', mouse, ' (output of python script)'])
    modalities = np.empty([0, 0])

    for i in range(numsess):
        sessdate = session_dates[i]
        print(sessdate)
        sessiondir = os.path.join(datadir,sessdate)
        filepath = glob.glob(os.path.join(sessiondir,'*_Session*.mat'))[0] #get filepath to behavioral data
        numtrials, numcorrect, stimtype, stimside, percentcorrect, correcttrials, choice, stimdict, targstim, diststim, singlespout, assisted, stimtypemean, rewarded, optotype = _getSessionData(filepath)
        writestring = ["\n\n", sessdate, " ------ session index ",str(i),
                       "\n\nAssisted = ",str(np.mean(assisted)),
                       "\nSinglespout = ",str(np.mean(singlespout)),
                       "\nDiststim mean (discrimination) = ",str(np.mean(diststim)),
                       "\nModality = ",str(np.mean(stimtypemean)),
                       "\nNumber of trials = ", str(numtrials),
                       "\nNumber of correct trials = ", str(numcorrect),
                       "\nPercent trials correct = ", str(percentcorrect)]
        file.writelines(writestring) #write the data above to text file
        modalities = np.append(modalities,stimtypemean)
    modalities = modalities.tolist()
    file.close()
    
    return session_dates, modalities

def _getSessionData(filepath):
    """
    Get data from a single session, add modality, assisted, singlespout, etc. 

    Parameters
    ----------
    filepath : TYPE
        Path to the desired behavioral (.mat) file for the session.

    Returns
    -------
    numtrials : TYPE
        DESCRIPTION.
    numcorrect : TYPE
        DESCRIPTION.
    stimtypes : TYPE
        DESCRIPTION.
    stimside : TYPE
        DESCRIPTION.
    percentcorrect : TYPE
        DESCRIPTION.
    correcttrials : TYPE
        DESCRIPTION.
    choice : TYPE
        DESCRIPTION.
    stimdict : TYPE
        DESCRIPTION.
    targstim : TYPE
        DESCRIPTION.
    diststim : TYPE
        DESCRIPTION.
    singlespout : TYPE
        DESCRIPTION.
    assisted : TYPE
        DESCRIPTION.
    stimtypemean : TYPE
        DESCRIPTION.
    rewarded : TYPE
        DESCRIPTION.
    optotype : TYPE
        DESCRIPTION.

    """
    importfile = sio.loadmat(filepath)
    choice = importfile['SessionData']['ResponseSide'][0, 0] - 1
    
    numtrials = int(importfile['SessionData']['nTrials'])
    stimside = importfile['SessionData']['CorrectSide'][0, 0] - 1
    stimtypes = importfile['SessionData']['StimType'][0, 0]
    stimstructs = importfile['SessionData']['stimEvents'][0,0]
    targstim = importfile['SessionData']['TargStim'][0,0]
    diststim = importfile['SessionData']['DistStim'][0,0]
    assisted = importfile['SessionData']['Assisted'][0,0]
    singlespout = importfile['SessionData']['SingleSpout'][0,0]
    rewarded = importfile['SessionData']['Rewarded'][0,0]
    try:
        optotype = importfile['SessionData']['optoType'][0,0] #get the opto data if it exists
    except ValueError:
        #print("No opto in this session")
        optotype = np.empty((1,numtrials))
        optotype[:,:] = np.nan #make opto data NaN if it doesn't exist
    
    correcttrials = np.multiply(1, choice == stimside)
    percentcorrect = round(np.sum(correcttrials)/numtrials, 3)
    numcorrect = np.sum(correcttrials)
    stimtypemean = np.mean(stimtypes)
    
    allleftclicks = []
    allrightclicks = []
    allleftflashes = []
    allrightflashes = []
    
    for i in range(numtrials): 
        stimevent = stimstructs[0,i]
        numleftclicks = np.size(stimevent[0,0])
        numrightclicks = np.size(stimevent[0,1])
        numleftflashes = np.size(stimevent[0,2])
        numrightflashes = np.size(stimevent[0,3])
        
        
        allleftclicks.append(numleftclicks)
        allrightclicks.append(numrightclicks)
        allleftflashes.append(numleftflashes)
        allrightflashes.append(numrightflashes)
        
    #print('\n' + filepath[39:] + " left click rate: " + str(rateleftclicks))
    #print('\n' + filepath[39:] + " right click rate: " + str(raterightclicks))
    
    stimdict = {'laudio': allleftclicks, 'raudio': allrightclicks, 'lvis': allleftflashes, 'rvis': allrightflashes}
    
    
    return numtrials, numcorrect, stimtypes, stimside, percentcorrect, correcttrials, choice, stimdict, targstim, diststim, singlespout, assisted, stimtypemean, rewarded, optotype

def formatSessions(filelist,input_dim,modality):
    """
    Calls _getSessionData for each session in filelist. Then puts that data into the necessary format for the 
    SSM package to train a GLM-Hmm

    Parameters
    ----------
    filelist : TYPE
        DESCRIPTION.
    input_dim : TYPE
        DESCRIPTION.
    modality : TYPE
        DESCRIPTION.

    Returns
    -------
    used_sessions : TYPE
        DESCRIPTION.
    allchoices : TYPE
        DESCRIPTION.
    allstimsides : TYPE
        DESCRIPTION.
    allnumtrials : TYPE
        DESCRIPTION.
    inpts : TYPE
        DESCRIPTION.
    alltargstim : TYPE
        DESCRIPTION.
    alldiststim : TYPE
        DESCRIPTION.

    """
    allchoices = []
    allstimsides = []
    allnumtrials = []
    inpts = []
    alltargstim = []
    alldiststim = []
    num_sess = len(filelist)
    used_sessions = np.ones(num_sess)
    print('\nRemoving trials with no choice on current trial and with opto stim\n\n')
    for i in range(num_sess):
        numtrials, numcorrect, stimtype, stimside, percentcorrect, correcttrials, choice, stimdict, targstim, diststim, singlespout, assisted, stimtypemean, rewarded, optotype = _getSessionData(filelist[i])
        
        if modality == 'audio':
            desired_modality = np.array(stimtype == 2)
        elif modality == 'visual':
            desired_modality = np.array(stimtype == 1)
        elif modality == 'tactile':
            print('Not working yet for tactile')
        elif modality == 'all':
            desired_modality = np.ones(numtrials) #array of ones, signifying every trial is in the "desired modality"
            
        if np.sum(desired_modality) / numtrials < .9: #if less than 90 percent of trials are in the desired modality
            used_sessions[i] = 0 #denote that this session is not used for training
            continue #skip this session in the for loop
        
        targstim = np.int32(targstim)
        diststim = np.int32(diststim)
        
        nochoicetrials = np.isnan(choice) # remove trials with no choice and with opto
        nooptotrials = np.isnan(optotype)
        choice = choice[np.logical_and(~nochoicetrials,nooptotrials)]
        choice = np.int32(choice)
        choice = np.expand_dims(choice, axis=1)
        
        stimside = stimside[np.logical_and(~nochoicetrials,nooptotrials)]
        targstim = targstim[np.logical_and(~nochoicetrials,nooptotrials)]
        diststim = diststim[np.logical_and(~nochoicetrials,nooptotrials)]
        
        lstim = np.zeros(len(stimside))
        rstim = np.zeros(len(stimside))
        lstim[stimside==0] = targstim[stimside==0]
        lstim[stimside==1] = diststim[stimside==1]
        rstim[stimside==0] = diststim[stimside==0]
        rstim[stimside==1] = targstim[stimside==1]
        
        numtrials_withchoice = choice.size  #number of trials per session with choice
        
        allchoices.append(choice)
        allstimsides.append(stimside)
        allnumtrials.append(numtrials_withchoice)
        alltargstim.append(targstim)
        alldiststim.append(diststim)
        
        temp = np.negative(np.ones([numtrials_withchoice,input_dim]))
        
        # #this is using stimdict, not targstim
        # #right now the following line is sign flipped
        # #coherence = np.divide(np.array(stimdict['laudio']), np.array(stimdict['laudio']) + np.array(stimdict['raudio']))#only audio trials right now
        # coherence = np.array(stimdict['laudio']) - np.array(stimdict['raudio'])
        # coherence = coherence[np.squeeze(np.logical_and(~nochoicetrials,nooptotrials))] #remove trials without choice and opto
        # #coherence = np.subtract(coherence,.5) #zero center
        # #coherence = coherence / np.std(coherence) #normalize
        # coherence = coherence / np.max(coherence) #normalize
        # temp[:,0] = coherence 
        
        #this is using targstim
        #right now the following line is sign flipped
        #coherence = np.divide(lstim, lstim + rstim)#only audio trials right now
        #coherence = np.subtract(coherence,.5) #zero center
        coherence = lstim - rstim
        #coherence = coherence / np.std(coherence) #normalize
        coherence = coherence / np.max(coherence) #normalize
        temp[:,0] = coherence 
        
        # #this separates out right and left
        # rstim = rstim / np.max(rstim)
        # lstim = lstim / np.max(lstim)
        # temp[:,0] = -rstim
        # temp[:,1] = lstim
        
        inpts.append(temp)
        
    print('From the ', str(num_sess), ' selected sessions, there are '+str(len(allchoices))+' that match the desired sensory modality and will be used to train the model\n')
    return used_sessions, allchoices, allstimsides, allnumtrials, inpts, alltargstim, alldiststim


if __name__ == '__main__':
    #dpath = 'X:\Widefield' #these lines must be commented out to work run this within MATLAB
    #mouse = 'CSP38'
    session_dates,modalities = createSessionTxtFile(dpath,mouse) #need to specify these args
