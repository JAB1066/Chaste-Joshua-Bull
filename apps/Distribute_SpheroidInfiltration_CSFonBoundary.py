import multiprocessing
import os
import subprocess
import itertools
import time

import numpy as np

executable = '/mi/share/scratch/bull/ChasteStuff/chaste-release/projects/JoshuaBull/apps/Exe_SpheroidInfiltration_CSFonBoundary'

chaste_test_dir = os.environ.get('CHASTE_TEST_OUTPUT')
path_to_output = os.path.join(chaste_test_dir,'ParameterSweeps','DorieReproduction','VaryCSF1ConcOnBoundary','Aug_8_2')

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

command_line_args = [' --ID ', ' --CC ', ' --ACCL ', ' --HC ', ' --QC ', ' --CHD ',' --AD ',' --OCR ',' --TTAM ',' --BCCSF ',' --IT ']
params_list = ['simulation_id', 'chemotaxisCoefficient', 'averageCellCycleLength', 'hypoxicConcentration', 'quiescentConcentration', 'criticalHypoxicDuration','apoptosisDuration','oxygenConsumptionRate','timeToAddMacrophages','boundaryConcCSF1', 'iterationNumber']

today = time.strftime('%Y-%m-%dT%H%M')

# Param ranges (in lists, for itertools product)
#accl = np.linspace(16, 32, num=2)
#hc = np.linspace(0.1, 0.7, num=4)
#qc = np.linspace(0.3, 0.7, num=3)
#chd = np.linspace(8.0, 16.0, num=2)
#ad = np.linspace(48.0, 48.0, num=1)
#ocr = np.linspace(0.03, 0.03, num=1)
#it = np.linspace(1, 10, num=10)
cc = np.linspace(10.0, 20.0, num=5)
accl = np.linspace(16, 16, num=1)
hc = np.linspace(0.3, 0.3, num=1)
qc = np.linspace(0.3, 0.3, num=1)
chd = np.linspace(8.0, 8.0, num=1)
ad = np.linspace(48.0, 48.0, num=1)
ocr = np.linspace(0.03, 0.03, num=1)
ttam = np.linspace(100, 300, num=5)
bccsf = np.linspace(0, 0, num=1)
it = np.linspace(1, 3, num=3)


combined_iterable = enumerate(itertools.product(cc,accl,hc,qc,chd,ad,ocr,ttam,bccsf,it))
shift = 1125 # For running large param sweeps, specify an offset if needed

def main():
    run_simulations()


# Create a list of commands and pass them to separate processes
def run_simulations():

    # Make a list of calls to a Chaste executable
    command_list = []
    
    if not os.path.exists(path_to_output):
    	os.makedirs(path_to_output)
    
    params_file = open(path_to_output + '/params_file.csv', 'w')
    print("Parameter file created at path: " + path_to_output + '/params_file.csv')
    params_file.write(','.join(params_list) + '\n')

    base_command = 'nice ' + executable

    for idx, param_set in combined_iterable:
    	idx = idx + shift
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
    # Wait at most one week (Now 10 weeks!)
    pool.map_async(execute_command, command_list).get(6048000)


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
