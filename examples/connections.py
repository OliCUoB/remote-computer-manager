from base_connection import Connection
import subprocess
import os
import re
import datetime
import pandas as pd
import subprocess
class PbsExample(Connection):
    """
    This is a mock example of the class one may write to utilise the remote-computer-management library. More specifically this inherits from the Connection class in base_connection.py.

    This is meant to conatin the BASIC commands that can be used by programs to control the remote computer.

    This example is to connect to a torque cluster with a PBS queuing system. The login shell will understand BASH.
    """

    def __init__(self, cluster_user_name, ssh_config_alias, path_to_key, forename_of_user, surname_of_user, user_email, base_output_path = '/base/output/path', base_runfiles_path = '/base/run/file/path', master_dir = '/master/dir'):
		Connection.__init__(self, cluster_user_name, ssh_config_alias, path_to_key, forename_of_user, surname_of_user, user_email)
		self.submit_command = 'qsub'
		self.information_about_cluster = 'Example Cluster Name (ECN): Advanced Computing Research Centre, somewhere.'
		self.base_output_path = base_output_path
		self.base_runfiles_path = base_runfiles_path
		self.master_dir = master_dir
		self.activate_venv_list = ['module add python-anaconda-4.2-3.5', 'source activate virtual_environment_name']

	# INSTANCE METHODS
	def checkQueue(self, job_number):
		"""This function must exist to satisfy the abstract class that it inherits from. In this case it takes a job number and returns a list of all the array numbers of that job still running."""
		grep_part_of_cmd = "qstat -tu " + self.user_name + " | grep \'" + str(job_number) + "\' | awk \'{print $1}\' | awk -F \"[][]\" \'{print $2}\'"

		output_dict = self.checkSuccess(self.sendCommand([grep_part_of_cmd])) # Remember that all commands should be passed through the "checkSuccess" function that is inherited from the Connection class.

		return output_dict

	def checkDiskUsage(self):
		"""This function returns disk usage details. In this case the cluster has a custom quota function 'quota_name'."""
		# create all the post connection commands needed
		get_disk_usage_units_command = "quota_name | awk \'{print $1}\' | tail -n 2 | head -n 1 | sed \'s/[<>]//g\'"
		get_disk_usage_command = "quota_name | awk \'{print $1}\' | tail -n 1"
		get_disk_usage_soft_limit_command = "quota_name | awk \'{print $2}\' | tail -n 1"
		get_disk_usage_hard_limit_command = "quota_name | awk \'{print $3}\' | tail -n 1"
		# combine the connection command with the post connection commands in a list (as is recomended).
		units_cmd = ["ssh", self.ssh_config_alias, get_disk_usage_units_command]
		usage_cmd = ["ssh", self.ssh_config_alias, get_disk_usage_command]
		soft_limit_cmd = ["ssh", self.ssh_config_alias, get_disk_usage_soft_limit_command]
		hard_limit_cmd = ["ssh", self.ssh_config_alias, get_disk_usage_hard_limit_command]
		# send the commands and save the exit codes and outputs
		units = Connection.getOutput(units_cmd)
		usage = Connection.getOutput(usage_cmd)
		soft_limit = Connection.getOutput(soft_limit_cmd)
		hard_limit = Connection.getOutput(hard_limit_cmd)
		# convert string outputs to floats where neccessary
		units[1] = str(units[1], "utf-8").rstrip()
		usage[1] = float(usage[1])
		soft_limit[1] = float(soft_limit[1])
		hard_limit[1] = float(hard_limit[1])
		# print some stats
		print(100 * (usage[1] / (1.0 * hard_limit[1]) ),"% of total disk space used.\n\n",hard_limit[1] - usage[1]," ",units[1]," left until hard limit.\n\n",soft_limit[1] - usage[1]," ",units[1]," left unit soft limit.", sep='')
		
		return usage, soft_limit, hard_limit, units

	def createStandardSubmissionScript(self, output_filename, pbs_job_name, queue_name, no_of_unique_kos, path_and_name_of_unique_ko_dir_names, no_of_repetitions_of_each_ko, wholecell_model_master_dir, output_dir, path_and_name_of_ko_codes, outfile_name_and_path, errorfile_name_and_path):
		# set job array numbers to None so that we can check stuff has wprked later
		job_array_numbers = None
		# The maximum job array size on BC3
		max_job_array_size = 500
		# initialise output dict
		output_dict = {}
		# test that a reasonable amount of jobs has been submitted (This is not a hard and fast rule but there has to be a max and my intuition suggestss that it will start to get complicated around this level i.e. queueing and harddisk space etc)
		total_sims = no_of_unique_kos * no_of_repetitions_of_each_ko
		if total_sims > 20000:
			raise ValueError('Total amount of simulations for one batch submission must be less than 20,000, here total_sims=',total_sims)

		output_dict['total_sims'] = total_sims
		# spread simulations across array jobs
		if no_of_unique_kos <= max_job_array_size:
			no_of_unique_kos_per_array_job = 1
			no_of_arrays = no_of_unique_kos
			job_array_numbers = '1-' + str(no_of_unique_kos)
			walltime = '30:00:00'
		else:
			# job_array_size * no_of_unique_kos_per_array_job = no_of_unique_kos so all the factors of no_of_unique_kos is
			common_factors = [x for x in range(1, no_of_unique_kos+1) if no_of_unique_kos % x == 0]
			# make the job_array_size as large as possible such that it is less than max_job_array_size
			factor_idx = len(common_factors) - 1
			while factor_idx >= 0:
				if common_factors[factor_idx] < max_job_array_size:
					job_array_numbers = '1-' + str(common_factors[factor_idx])
					no_of_arrays = common_factors[factor_idx]
					no_of_unique_kos_per_array_job = common_factors[(len(common_factors)-1) - factor_idx]
					factor_idx = -1
				else:
					factor_idx -= 1

			# raise error if no suitable factors found!
			if job_array_numbers is None:
				raise ValueError('job_array_numbers should have been assigned by now! This suggests that it wasn\'t possible for my algorithm to split the KOs across the job array properly. Here no_of_unique_kos=', no_of_unique_kos, ' and the common factors of this number are:', common_factors)

			# add some time to the walltime because I don't think the jobs have to startat the same time
			walltime = '35:00:00'

		output_dict['no_of_arrays'] = no_of_arrays
		output_dict['no_of_unique_kos_per_array_job'] = no_of_unique_kos_per_array_job
		output_dict['no_of_repetitions_of_each_ko'] = no_of_repetitions_of_each_ko
		# calculate the amount of cores per array job - NOTE: for simplification we only use cores and not nodes (this is generally the fastest way to get through the queue anyway)
		no_of_cores = no_of_repetitions_of_each_ko * no_of_unique_kos_per_array_job
		output_dict['no_of_sims_per_array_job'] = no_of_cores
		output_dict['list_of_rep_dir_names'] = list(range(1, no_of_repetitions_of_each_ko + 1))
		no_of_nodes = 1
		# write the script to file
		with open(output_filename, mode='wt', encoding='utf-8') as myfile:
			myfile.write("#!/bin/bash" + "\n")
			myfile.write("\n")
			myfile.write("# This script was automatically created by Oliver Chalkley's whole-cell modelling suite. Please contact on o.chalkley@bristol.ac.uk" + "\n")
			myfile.write("# Title: " + pbs_job_name + "\n")
			myfile.write("# User: " + self.forename_of_user + ", " + self.surename_of_user + ", " + self.user_email + "\n")
			myfile.write("# Affiliation: Minimal Genome Group, Life Sciences, University of Bristol " + "\n")
			myfile.write("# Last Updated: " + str(datetime.datetime.now()) + "\n")
			myfile.write("\n")
			myfile.write("# BC3: 223 base blades which have 16 x 2.6 GHz SandyBridge cores, 4GB/core and a 1TB SATA disk." + "\n")
			myfile.write("\n")
			myfile.write("## Job name" + "\n")
			myfile.write("#PBS -N " + pbs_job_name + "\n")
			myfile.write("\n")
			myfile.write("## Resource request" + "\n")
			myfile.write("#PBS -l nodes=" + str(no_of_nodes) + ":ppn=" + str(no_of_cores) + ",walltime=" + walltime + "\n")
			myfile.write("#PBS -q " + queue_name + "\n")
			myfile.write("\n")
			myfile.write("## Job array request" + "\n")
			myfile.write("#PBS -t " + job_array_numbers + "\n")
			myfile.write("\n")
			myfile.write("## designate output and error files" + "\n")
			myfile.write("#PBS -e " + outfile_name_and_path + "\n")
			myfile.write("#PBS -o " + errorfile_name_and_path + "\n")
			myfile.write("\n")
			myfile.write("# print some details about the job" + "\n")
			myfile.write('echo "The Array ID is: ${PBS_ARRAYID}"' + "\n")
			myfile.write('echo Running on host `hostname`' + "\n")
			myfile.write('echo Time is `date`' + "\n")
			myfile.write('echo Directory is `pwd`' + "\n")
			myfile.write('echo PBS job ID is ${PBS_JOBID}' + "\n")
			myfile.write('echo This job runs on the following nodes:' + "\n")
			myfile.write('echo `cat $PBS_NODEFILE | uniq`' + "\n")
			myfile.write("\n")
			myfile.write("# load required modules" + "\n")
			myfile.write("module unload apps/matlab-r2013b" + "\n")
			myfile.write("module load apps/matlab-r2013a" + "\n")
			myfile.write('echo "Modules loaded:"' + "\n")
			myfile.write("module list" + "\n")
			myfile.write("\n")
			myfile.write("# create the master directory variable" + "\n")
			myfile.write("master=" + wholecell_model_master_dir + "\n")
			myfile.write("\n")
			myfile.write("# create output directory" + "\n")
			myfile.write("base_outDir=" + output_dir + "\n")
			myfile.write("\n")
			myfile.write("# collect the KO combos" + "\n")
			myfile.write("ko_list=" + path_and_name_of_ko_codes + "\n")
			myfile.write("ko_dir_names=" + path_and_name_of_unique_ko_dir_names + "\n")
			myfile.write("\n")
			myfile.write("# Get all the gene KOs and output folder names" + "\n")
			myfile.write('for i in `seq 1 ' + str(no_of_unique_kos_per_array_job) + '`' + "\n")
			myfile.write('do' + "\n")
			myfile.write('    Gene[${i}]=$(awk NR==$((' + str(no_of_unique_kos_per_array_job) + '*(${PBS_ARRAYID}-1)+${i})) ${ko_list})' + "\n")
			myfile.write('    unique_ko_dir_name[${i}]=$(awk NR==$((' + str(no_of_unique_kos_per_array_job) + '*(${PBS_ARRAYID}-1)+${i})) ${ko_dir_names})' + "\n")
			myfile.write("done" + "\n")
			myfile.write("\n")
			myfile.write("# go to master directory" + "\n")
			myfile.write("cd ${master}" + "\n")
			myfile.write("\n")
			myfile.write("# NB have limited MATLAB to a single thread" + "\n")
			myfile.write('options="-nodesktop -noFigureWindows -nosplash -singleCompThread"' + "\n")
			myfile.write("\n")
			myfile.write("# run 16 simulations in parallel")
			myfile.write('echo "Running simulations (single threaded) in parallel - let\'s start the timer!"' + "\n")
			myfile.write('start=`date +%s`' + "\n")
			myfile.write("\n")
			myfile.write("# create all the directories for the diarys (the normal output will be all mixed up cause it's in parrallel!)" + "\n")
			myfile.write('for i in `seq 1 ' + str(no_of_unique_kos_per_array_job) + '`' + "\n")
			myfile.write("do" + "\n")
			myfile.write('    for j in `seq 1 ' + str(no_of_repetitions_of_each_ko) + '`' + "\n")
			myfile.write("    do" + "\n")
			myfile.write('        specific_ko="$(echo ${Gene[${i}]} | sed \'s/{//g\' | sed \'s/}//g\' | sed \"s/\'//g\" | sed \'s/\"//g\' | sed \'s/,/-/g\')/${j}"' + "\n")
			myfile.write('        mkdir -p ${base_outDir}/${unique_ko_dir_name[${i}]}/diary${j}' + "\n")
			myfile.write('        matlab ${options} -r "diary(\'${base_outDir}/${unique_ko_dir_name[${i}]}/diary${j}/diary.out\');addpath(\'${master}\');setWarnings();setPath();runSimulation(\'runner\',\'koRunner\',\'logToDisk\',true,\'outDir\',\'${base_outDir}/${unique_ko_dir_name[${i}]}/${j}\',\'jobNumber\',$((no_of_repetitions_of_each_ko*no_of_unique_kos_per_array_job*(${PBS_ARRAYID}-1)+no_of_unique_kos_per_array_job*(${i}-1)+${j})),\'koList\',{{${Gene[${i}]}}});diary off;exit;" &' + "\n")
			myfile.write("    done" + "\n")
			myfile.write("done" + "\n")
			myfile.write("wait" + "\n")
			myfile.write("\n")
			myfile.write("end=`date +%s`" + "\n")
			myfile.write("runtime=$((end-start))" + "\n")
			myfile.write('echo "$((${no_of_unique_kos_per_array_job}*${no_of_repetitions_of_each_ko})) simulations took: ${runtime} seconds."')

		# give the file execute permissions
		subprocess.check_call(["chmod", "700", str(output_filename)])

		return output_dict

	def getJobIdFromSubStdOut(self, stdout):
		return int(re.search(r'\d+', stdout).group())
