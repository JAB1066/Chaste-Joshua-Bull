import multiprocessing
import os
import subprocess
import itertools
import time

import numpy as np

executable = '/mi/share/scratch/bull/ChasteStuff/chaste-release/projects/JoshuaBull/apps/Exe_ParallelDorieReproduction_SteadyStateQuest'

chaste_test_dir = os.environ.get('CHASTE_TEST_OUTPUT')
path_to_output = os.path.join(chaste_test_dir,'ParameterSweeps','DorieReproduction_SteadyStateQuest_2D')

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

command_line_args = [' --ID ', ' --ACCL ', ' --HC ', ' --QC ', ' --CHD ']
params_list = ['simulation_id', 'averageCellCycleLength', 'hypoxicConcentration', 'quiescentConcentration', 'criticalHypoxicDuration']

today = time.strftime('%Y-%m-%dT%H%M')

# Param ranges (in lists, for itertools product)
accl = np.linspace(8.0, 14.0, num=4)
hc = np.linspace(0.05, 0.35, num=7)
qc = np.linspace(0.35, 0.6, num=6)
chd = np.linspace(2.0, 10.0, num=5)

combined_iterable = enumerate(itertools.product(accl,hc,qc,chd))


def main():
    run_simulations()


# Create a list of commands and pass them to separate processes
def run_simulations():

    # Make a list of calls to a Chaste executable
    command_list = []
    
    if not os.path.exists(path_to_output):
    	os.makedirs(path_to_output)
    
    params_file = open(path_to_output + '/params_file.csv', 'w')
    params_file.write(','.join(params_list) + '\n')

    base_command = 'nice ' + executable

    for idx, param_set in combined_iterable:
    
    	params_file.write(str(idx) + ',' + ",".join(map(str, param_set)) + '\n')

        command = base_command 
        command += ' --ID ' + str(idx)
        
        for arg in range(len(param_set)):
        	command += command_line_args[arg+1] + str(param_set[arg])

        command_list.append(command)
        
    params_file.close()

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    # Wait at most two weeks
    pool.map_async(execute_command, command_list).get(1209600)


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
