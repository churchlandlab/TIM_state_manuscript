import numpy as np
import os
import sys


def create_params_array(savepath, animals, methods, N_iters, N_states, alphas, sigmas): #must be list or numpy array
    if os.path.isfile(os.path.join(savepath, 'params.npy')):
        print("Params file already exist in this folder. Skipping creation.")
        param_matrix = np.load(os.path.join(savepath,'params.npy'), allow_pickle=True)
        return param_matrix.shape[0], param_matrix.shape[1], param_matrix

    NUM_PARAMS = 6
    param_matrix = np.empty((0,NUM_PARAMS))

    for animal in animals:
        for mthd in methods:
            for N_iter in N_iters:
                for N_state in N_states:
                    for alpha in alphas:
                        for sigma in sigmas:
                            paramset = np.array([[animal, mthd, N_iter, N_state, alpha, sigma]], dtype=object) # must be 2D for np.append
                            param_matrix = np.append(param_matrix, paramset, axis=0)
    
    numruns = param_matrix.shape[0]
    numparams = param_matrix.shape[1]
    print(numruns) #send to stdout for bash variable  

    if not os.path.exists(savepath):
        os.makedirs(savepath)
        
    paramfile = os.path.join(savepath, 'params.npy')
    with open(paramfile, 'wb') as f:
        np.save(f, param_matrix)
        
    jobnumfile = os.path.join(savepath, 'jobnum.npy')
    with open(jobnumfile, 'wb') as f:
        np.save(f, numruns)

    return numruns, numparams, param_matrix

def build_submission_script(savepath, filename, pyfile, numjobs, n_cpu_per_job=3, pyfile_args=[]): #pyfile is the absolute path of the python file to be run
    # args should be a list of args that will be passed to the python file we wish to run
    if os.path.isfile(os.path.join(savepath, filename)):
        raise FileExistsError("Submission script already exists")
    
    submission_text = '''#!/bin/bash

mkdir -p {savepath}

#$ -cwd
# error = Merged with joblog
#$ -o {savepath}/logs/joblog.$JOB_ID.$TASK_ID

#$ -j y
#$ -l h_rt=24:00:00,h_data=4G
#$ -pe shared {n_cpu_per_job}

#$ -t 1-{numjobs}:1

. /u/local/Modules/default/init/modules.sh
source $HOME/.bash_profile
module load anaconda3
conda activate max_glmhmm


echo "JOB started on:  " `hostname -s`
echo "JOB started on:  " `date `
echo " "
echo " " # the command goes in quotes   


echo "This is sub-job $SGE_TASK_ID"
echo "{pyfile} {args} "$SGE_TASK_ID""

python {pyfile} {args} "$SGE_TASK_ID" #run the python script with the proper index

echo "JOB ended on:  " `hostname -s`
echo "JOB ended on:  " `date `
echo " "
'''.format(savepath=savepath,
           numjobs=numjobs,
           pyfile=pyfile,
           n_cpu_per_job=n_cpu_per_job,
           args=' '.join(str(i) for i in pyfile_args))

    fname = os.path.join(savepath, f'{filename}') 
    with open(fname,'w') as file:
        file.writelines(submission_text)
    return fname


def stitch_outputs(savepath):
    
    jobnum = np.load(os.path.join(savepath,'jobnum.npy'))
    searchdir = os.path.join(savepath,'outputs')
    files = sorted(os.listdir(searchdir))
    numoutputs = len(np.load(os.path.join(searchdir, files[1]))) #load one file to determine how many outputs
    print('There are {} output files found. There are {} output files expected.'.format(len(files), jobnum))

    #outputs = np.empty((jobnum,numoutputs))
    #outputs[:] = np.nan
    outputs = [None] * int(jobnum)
    for file in files:
        out = np.load(os.path.join(searchdir, file))
        index = int(file.replace('.npy',''))
        #outputs[index,:] = out
        outputs[index] = out
    return outputs


if __name__ == '__main__':
    # savepath = '/u/home/m/mmelin/UGE_GLM_runs/run1' #for debugging
    # stitch_outputs(savepath,864)
    
    savepath = sys.argv[1] #path for this particular run
    numruns = create_params_array(savepath)
