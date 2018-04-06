import multiprocessing
import os
import subprocess
import itertools
import time

import numpy as np

executable = '/mi/share/scratch/bull/ChasteStuff/chaste-release/projects/JoshuaBull/apps/Exe_ParallelFullSpheroidSimulationsWithMacrophages'
#export CHASTE_TEST_OUTPUT=/media/bull/MyPassport/DPhil/chaste_test_output/
#export CHASTE_TEST_OUTPUT=/mi/share/scratch/bull/ChasteStuff/chaste_test_output/
chaste_test_dir = os.environ.get('CHASTE_TEST_OUTPUT')

#path_to_output = os.path.join(chaste_test_dir, 'FullSimulations','DiffusionOnlyParameterSweep')
path_to_output = os.path.join(chaste_test_dir, 'ParameterSweeps','Dorie1982Reproduction','VaryBrownianMotion')

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

command_line_args = [' --ID ', ' --DC ', ' --IN ']
params_list = ['simulation_id', 'diffusionCoefficient','iterationNumber']

today = time.strftime('%Y-%m-%dT%H%M')

# Param ranges (in lists, for itertools product)
dc = np.linspace(0.0, 2.0, num=41)
iterationNumber = np.linspace(1.0,25.0,num=25)

combined_iterable = enumerate(itertools.product(dc,iterationNumber))


def main():
    run_simulations()


# Create a list of commands and pass them to separate processes
def run_simulations():

    # Make a list of calls to a Chaste executable
    command_list = []
    
    if not os.path.exists(path_to_output):
    	os.makedirs(path_to_output)
    
    params_file = open(path_to_output + '/params_file_wolverine.csv', 'w')
    params_file.write(','.join(params_list) + '\n')

    base_command = executable

    for idx, param_set in combined_iterable:
    
    	params_file.write(str(idx) + ',' + ",".join(map(str, param_set)) + '\n')

        command = base_command 
        command += ' --ID ' + str(idx)
        
        for arg in range(len(param_set)):
        	command += command_line_args[arg+1] + str(param_set[arg])

        command_list.append(command)
        
    params_file.close()

    # Use processes equal to the number of cpus available
    count = 24#multiprocessing.cpu_count()

    print("Py: Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    # Wait at most one day
    pool.map_async(execute_command, command_list).get(86400)


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
