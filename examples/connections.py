from base_connection import Connection
import subprocess
import re
import datetime
class Pbs(Connection):
    """
    This is meant to be a template to create a connection object for a standard PBS/TORQUE cluster. This inherits from the Connection class in base_connection.py.

    This is meant to contain the BASIC commands that can be used by programs to control the remote computer (that aren't already included in base_connection.Connection). This is atomistic level commands that form the basis of more complex and specific programs.
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
        """
        This function must exist to satisfy the abstract class that it inherits from. In this case it takes a job number and returns a list of all the array numbers of that job still running.
        
        Args:
            job_number (int): PBS assigns a unique integer number to each job. Remeber that a job can actually be an array of jobs.

        Returns:
            output_dict (dict): Has keys 'return_code', 'stdout', and 'stderr'.
        """

                # -t flag shows all array jobs related to one job number, if that job is an array.
        grep_part_of_cmd = "qstat -tu " + self.user_name + " | grep \'" + str(job_number) + "\' | awk \'{print $1}\' | awk -F \"[][]\" \'{print $2}\'"

        output_dict = self.checkSuccess(self.sendCommand([grep_part_of_cmd])) # Remember that all commands should be passed through the "checkSuccess" function that is inherited from the Connection class.

        return output_dict

    def checkDiskUsage(self):
        """
        This function returns disk usage details of the remote computer. In this case the cluster has a custom quota function 'quota_name'.
        
        Returns:
            usage (float): The amount of disk space currently used.
            soft_limit (float): The soft limit (warnings if you go over) of the user's dis space.
            hard_limit (float): The hard limit (cannot write to disk at this point)
            units (str): The units of the data i.e. GB or TB etc.
        """

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

    def createStandardSubmissionScript(self, output_filename, pbs_job_name, queue_name, no_of_unique_jobs, no_of_repetitions_of_each_job, master_dir, outfile_name_and_path, errorfile_name_and_path, walltime, initial_message_in_code, list_of_job_specific_code):
        """
        This acts as a template for a submission script for the cluster however it does not contain any code for specific jobs. This code is pass to the function through the list_of_job_specific_code variable.

        The format for a submission in this case will be an array of jobs. Here we want to be able to specify a number of unique jobs and then the amount of times we wish to repeat each unique job. This will then split all the jobs across arrays and CPUs on the cluster depending on how many are given. Each unique job has a name and some settings, this is stored on the cluster in 2 files job_names_file and job_settings_file, respectively.

        Args:
            output_filename (str): The name of the submission script.
            pbs_job_name (str): The name given to the queuing system.
            queue_name (str): This cluster has a choice of queues and this variable specifies which one to use.
            no_of_unique_jobs (int): Total amount of jobs to run.
            no_of_repetitions_of_each_job (int): Total amount of repetitions of each job.
            master_dir (str): The directory on the remote computer that you want the submission script to start in.
            outfile_name_and_path (str): Absolute path and file name of where you want the outfiles of each job array stored.
            errorfile_name_and_path (str): Absolute path and file name of where you want to store the errorfiles of each job array stored.
            walltime (str): The maximum amount of time the job is allowed to take. Has the form 'HH:MM:SS'.
            initial_message_in_code (str): The first comment in the code normally says a little something about where this script came from. NOTE: You do not need to include a '#' to indicat it is a comment.
            list_of_job_specific_code (list of strings): Each element of the list contains a string of one line of code.

        Returns:
            output_dict (dict): Contains details of how it spread the jobs across arrays and CPUs. Has keys, 'no_of_arrays', 'no_of_unique_jobs_per_array_job', 'no_of_repetitions_of_each_job', 'no_of_sims_per_array_job', and 'list_of_rep_dir_names'.
        """

        # set job array numbers to None so that we can check stuff has worked later
        job_array_numbers = None
        # The maximum job array size on the cluster.
        max_job_array_size = 500
        # initialise output dict
        output_dict = {}
        # test that a reasonable amount of jobs has been submitted (This is not a hard and fast rule but there has to be a max and my intuition suggestss that it will start to get complicated around this level i.e. queueing and harddisk space etc)
        total_sims = no_of_unique_jobs * no_of_repetitions_of_each_job
        if total_sims > 20000:
            raise ValueError('Total amount of simulations for one batch submission must be less than 20,000, here total_sims=',total_sims)

        output_dict['total_sims'] = total_sims
        # spread simulations across array jobs
        if no_of_unique_jobs <= max_job_array_size:
            no_of_unique_jobs_per_array_job = 1
            no_of_arrays = no_of_unique_jobs
            job_array_numbers = '1-' + str(no_of_unique_jobs)
        else:
            # job_array_size * no_of_unique_jobs_per_array_job = no_of_unique_jobs so all the factors of no_of_unique_jobs is
            common_factors = [x for x in range(1, no_of_unique_jobs+1) if no_of_unique_jobs % x == 0]
            # make the job_array_size as large as possible such that it is less than max_job_array_size
            factor_idx = len(common_factors) - 1
            while factor_idx >= 0:
                if common_factors[factor_idx] < max_job_array_size:
                    job_array_numbers = '1-' + str(common_factors[factor_idx])
                    no_of_arrays = common_factors[factor_idx]
                    no_of_unique_jobs_per_array_job = common_factors[(len(common_factors)-1) - factor_idx]
                    factor_idx = -1
                else:
                    factor_idx -= 1

            # raise error if no suitable factors found!
            if job_array_numbers is None:
                raise ValueError('job_array_numbers should have been assigned by now! This suggests that it wasn\'t possible for my algorithm to split the KOs across the job array properly. Here no_of_unique_jobs=', no_of_unique_jobs, ' and the common factors of this number are:', common_factors)

        output_dict['no_of_arrays'] = no_of_arrays
        output_dict['no_of_unique_jobs_per_array_job'] = no_of_unique_jobs_per_array_job
        output_dict['no_of_repetitions_of_each_job'] = no_of_repetitions_of_each_job
        # calculate the amount of cores per array job - NOTE: for simplification we only use cores and not nodes (this is generally the fastest way to get through the queue anyway)
        no_of_cores = no_of_repetitions_of_each_job * no_of_unique_jobs_per_array_job
        output_dict['no_of_sims_per_array_job'] = no_of_cores
        output_dict['list_of_rep_dir_names'] = list(range(1, no_of_repetitions_of_each_job + 1))
        no_of_nodes = 1
        # write the script to file
        with open(output_filename, mode='wt', encoding='utf-8') as myfile:
            myfile.write("#!/bin/bash" + "\n")
            myfile.write("\n")
            myfile.write("# This script was created using Oliver Chalkley's remote-computer-manager library - https://github.com/OliCUoB/remote-computer-manager." + "\n")
            myfile.write("# " + initial_message_in_code + "\n")
            myfile.write("# Title: " + pbs_job_name + "\n")
            myfile.write("# User: " + self.forename_of_user + ", " + self.surename_of_user + ", " + self.user_email + "\n")
            if type(self.affiliation) is not None:
                myfile.write("# Affiliation: " + self.affiliation + "\n")
            myfile.write("# Last Updated: " + str(datetime.datetime.now()) + "\n")
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
            for line in list_of_job_specific_code:
                myfile.write(line)

        # give the file execute permissions
        subprocess.check_call(["chmod", "700", str(output_filename)])

        return output_dict

    def getJobIdFromSubStdOut(self, stdout):
        """
        When one submits a job to the cluster it returns the job ID to the stdout. This function takes that stdout and extracts the job ID so that it can be used to monitor the job if neccessary.

        Args:
            stdout (str): The stdout after submitting a job to the queue.

        Returns:
            return (int): The job ID of the job submitted which returned stdout.
        """
        
        return int(re.search(r'\d+', stdout).group())

    def createLocalFile(self):
        # createLocalFile needs to exist in order to satisfy the abstract method in the parent class. However in this example we decide that the only file we will ever need to create are job submission scripts and so simply pass the function.
        pass
