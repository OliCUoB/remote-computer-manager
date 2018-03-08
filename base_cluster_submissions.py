from abc import ABCMeta, abstractmethod
class BaseJobSubmission(metaclass=ABCMeta):
    """
    This is an abstract class that all job submission class should inherit from and acts as a template so that all higher level programs can take advantage of this structure and know that future clusters/submission queues will still work.

    This class assumes that the cluster connection takes the form of the base_connection class.
    """
    def __init__(self, submission_name, cluster_connection, simulation_output_path, errorfile_path, outfile_path, runfiles_path, repeitions_of_unique_task, master_dir, temp_storage_path):
        """
        The general idea of the structure is that all job submissions will require atleast a computer cluster, a job submission script (and all the details needed to make that script) and a command to submit the job to the cluster queuing system. This is meant to be as abstract/general as possible so things that are specific to a specific cluster should be included in a child class.

        Args:
            submission_name (str): A name that describes this job submission that can be used in job names and file names etc.
            cluster_connection (connection): This is a connection object that has inherited from the base_connection.Connection class.
            simulation_output_path (str): The absolute path on the cluster that you want any output from your program to be stored.
            errorfile_path (str): The absolute path on the cluster that you want the queuing systems error files to be stored.
            outfile_path (str): The absolute path on the cluster that you want the queuing systems output files to be stored.
            runfiles_path (str): The absolute path on the cluster where you want the job submission file stored.
            repeitions_of_unique_task (int): The number of times you want your program executed. Note: this is exact repetitions.
            master_dir (str): The absolute path on the cluster that you want the submission script to cd into.
            temp_storage_path (str): The absolute path on the local computer that you want temporary files to be stored on.
        """
        
        self.submission_name = submission_name
        self.submission_file_name = self.submission_name + '_sub_file.sh'
        self.temp_storage_path = temp_storage_path
        self.unique_job_name = self.createUniqueJobName(self.submission_name + '_')
        os.makedirs(self.temp_storage_path + '/' + self.unique_job_name)
        self.temp_storage_path = self.temp_storage_path + '/' + self.unique_job_name
        self.cluster_connection = cluster_connection
        self.simulation_output_path = simulation_output_path
        self.errorfile_path = errorfile_path
        self.outfile_path = outfile_path
        self.runfiles_path = runfiles_path
        self.submission_script_tmp_storage = self.temp_storage_path + '/' + self.submission_file_name
        self.list_of_directories_to_make_on_cluster = self.createListOfClusterDirectoriesNeeded()
        self.file_source_to_file_dest_dict = self.createDictOfFileSourceToFileDestinations()
        self.cluster_job_number = None
        self.time_of_submission = None

    # ABSTRACT METHODS
    ## createAllFiles function creates all the files needed by the submission. This will vary depending on type of job and so is left as an abstract method.
    ## The createListOfClusterDirectoriesNeeded and createDictOfFileSourceToFileDestinations functions are needed in order to know what files and where they need to be on the cluster. This will vary depending on type of job and so is left as an abstract method that will be properly codeed in the child class that will inherit from this.
    @abstractmethod
    def createAllFiles(self):
        pass

    @abstractmethod
    def createListOfClusterDirectoriesNeeded(self):
        pass

    @abstractmethod
    def createDictOfFileSourceToFileDestinations(self):
        pass

    def prepareForSubmission(self):
        """
        Makes sure that:
            1. All files needed are created.
            2. Any directories that will be needed by the job are present on the cluster.
            3. Transfers all neccessary files to the cluster.

        This information is taken from the relavent class variables which should be updated by the 'createAllFiles', 'createListOfClusterDirectoriesNeeded', and 'createDictOfFileSourceToFileDestinations' functions.

        Returns:
            return (list): This is a list of all the output dicts returrned from cluster connections. It is made up of the following:
                                - makedir_output_dict (dict): The output dict returned from the connection command to create the directories on the cluster.
                                - list_of_transferFiles_output_dicts (list of dicts): Each dict is the output dict returned from the connection command that transfered the files to the cluster. All transfers will have their output dict appended to a list and the list is returned.

        """

        # 1. All files needed are created.
        self.createAllFiles()

        # 2. construct the bash command to create the neccessary directories should they not be present.
        self.list_of_directories_to_make_on_cluster
        makedir_commands_list = ['mkdir -p ' + directory for directory in self.list_of_directories_to_make_on_cluster]
        # create directories on the cluster passing the "sendCommand" function through the "checkSuccess" function which are from the cluster_connection (this is of the form of the base_connection abstract class) which was passed when this class was instantiated.
        makedir_output_dict = self.cluster_connection.checkSuccess(self.cluster_connection.sendCommand, makedir_commands_list)

        # 3. transfer all neccessary files to the cluster using the appropriate functions from the cluster_connection instance.
        list_of_transferFiles_output_dicts = []
        for file_source in self.file_source_to_file_dest_dict.keys():
            list_of_transferFiles_output_dicts.append(self.cluster_connection.checkSuccess(self.cluster_connection.transferFile, file_source, self.file_source_to_file_dest_dict[file_source]))

        return [makedir_output_dict] + list_of_transferFiles_output_dicts

    def submitJobToCluster(self):
        """
        This function submits a job to the cluster queue, records the time that the connection returns it's output dict, retrieves the corresponding job number, and deletes all the local files created to make his submission happen.

        Returns:
            submit_job_ouput_dict (dict): The connection output dict returned once the submission was successfully executed.
        """

        # Create the job submission command
        submit_command = self.cluster_connection.submit_command + ' ' + self.runfiles_path + '/' + self.submission_file_name
        # Submit the job to the cluster queue
        submit_job_ouput_dict = self.cluster_connection.checkSuccess(self.cluster_connection.sendCommand, [submit_command])
        # Record the time that the connection returned it's output dict
        now = datetime.datetime.now()
        self.time_of_submission = {'day': now.day, 'month': now.month, 'year': now.year}
        # Record the job number of the submitted job
        self.cluster_job_number = self.cluster_connection.getJobIdFromSubStdOut(submit_job_ouput_dict['stdout'])
        # tidy up the tmp storage area if the submission was successful
        if submit_job_ouput_dict['return_code'] == 0:
            tmp_return_code = subprocess.call("rm -r " + self.temp_storage_path, shell=True)
            if tmp_return_code != 0:
                print("WARNING!!!! Could not remove the temporary files from ", self.temp_storage_path, " please fix this problem ASAP since if it carries on repeating it will filll the entire computer up until it breaks!")

        return submit_job_ouput_dict

    def createUniqueJobName(self, prefix):
        """
        Temporary files need a place to be stored. This creates a unique name that can be used to name the temporary directory. The name takes the form string + digits where the digits are created based on the current time.

        Args:
            prefix (str): The beginning part of the temporary directory name.

        Returns:
            return (str): prefix + unique_digits
        """
        valid_name = False
        while valid_name == False:
            unique_name_end = str('{:.9f}'.format(time.time()))
            unique_name_end = unique_name_end.replace(".","")
            unique_path = self.temp_storage_path + '/' + prefix + unique_name_end
            if os.path.isdir(unique_path) == False:
                valid_name = True

        return prefix + unique_name_end

class BaseManageSubmission(metaclass=ABCMeta):
    """
    This is an abstract class that all manage submission classes should inherit from so that it maintains a structure that higher level programs can reply on.

    The class initialises with a submission instance which will be a class that inherits from the BaseJobSubmission class and prepares the job for submission, submits the job to the cluster and then executes the abstract method 'monitorSubmission' to monitor the submission so that upper level programs can do things like waiting until a job is fnished or automatically do data processsing or updating databases with resultss etc.
    """
    def __init__(self, submission_instance):
        """
        This function takes a submission instance that is a class that inherits from the BaseJobSubmission class and prepares the job for submission, submits the job to the cluster and then monitors the jobs progress.

        Args:
            submission_instance (JobSubmission): The JobSubmission object is an object that inherits from the BaseJobSubmission class. 
        """
        self.submission = submission_instance
        outs =  self.submission.prepareForSubmission()
        list_of_return_codes = [output['return_code'] for output in outs]
        if sum(list_of_return_codes) > 0:
            raise ValueError('There has been a problem preparing some jobs for submission. The return codes from JobSubmission.prepareForSubmission() are: ', list_of_return_codes, '. The submission file name is ', self.submission.submission_file_name, ' and the temp_storage_path was ', self.submission.temp_storage_path)

        submission_out_dict = self.submission.submitJobToCluster()
        if submission_out_dict['return_code'] > 0:
            raise ValueError('There has been a problem submitting a job. The return code from JobSubmission.submitJobToCluster() is: ', submission_out_dict['return_code'], ". The submission file name is ", self.submission.submission_file_name, " and the temp_storage_path was ", self.submission.temp_storage_path)

        self.monitorSubmission(self.submission)

        return 

    # ABSTRACT METHODS
    @abstractmethod
    # This method is to monitor the progress of a job and perform other job related to the job like data processing and updating of databases etc
    def monitorSubmission(self, submission):
        pass
