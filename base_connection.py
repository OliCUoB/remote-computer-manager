from abc import ABCMeta, abstractmethod
import subprocess
import sys
import os
import time
class Connection(metaclass=ABCMeta):
    """
    This is an abstract class that all connection classes inherit from. The purpose of this class is to act as a template with which to communicate with other computers in a rigid manner so that other programs can be built on top of it, without knowing what computers it might connect to iin the future.
    
    REMEMBER that all connections to remote computers should go through the 'checkSuccess' function so that connection errors can be dealt with robustly.

    The general idea is that when communicating between computers it is often useful to be able to control them, send files to them and check how much storage space you need. If connecting to computing clusters you may also wish to be able to check the job queue.

    MOST IMPORTANTLY whenever you make any kind of connection to the remote computer you want to make sure that the remote computer got the messge and if it didn't then you need to be able to keep retrying without overloading the remote computer with connection requests (i.e. DDoS attack).
    """
    
    def __init__(self, remote_user_name, ssh_config_alias, path_to_key, forename_of_user, surname_of_user, user_email, affiliation = None):
        """
        In order to initiate this class the user must have their ssh config file set up to have their cluster connection as an alias. It is best to set this up on a secure ccomputer that you trusst and have an encryption key without a password. Details about setting up the SSH config file can be found at my website.
        
        To better explain things I will describe a toy example that all doc strings in this class will refer to. The user has set up an easy connecttion to the remote computer by setting up their ~/.ssh/config file (either directly or through a tunnel) like:

        Host ssh_alias
            User user_name
            HostName address_to_remote_computer
            IdentityFile /home/local_user_name/.ssh/path_to_key/key_name

        Args:
            remote_user_name (str): The username used on the remote computer, (i.e. user_name).
            ssh_config_alias (str): The name given to the SSH connection (i.e. ssh_config_alias).
            path_to_key (str): The path, and name, to the encryption key (i.e. /home/local_user_name/.ssh/path_to_key/key_name).
            forename_of_user (str): Your first name.
            surname_of_user (str): Your surname.
            user_email (str): Your email address.
        """

        self.user_name = remote_user_name
        self.ssh_config_alias = ssh_config_alias
        self.path_to_key = path_to_key
        self.forename_of_user = forename_of_user
        self.surename_of_user = surname_of_user
        self.user_email = user_email
        self.affiliation = affiliation

    # ABSTRACT METHODS
    @abstractmethod
    def createLocalFile(self):
        # Often one will want to create files on the local computer and transfer them to the remote computer. What files and how you create them depends on the specific task and so this is left as an abstract method.
        pass

    @abstractmethod
    def createStandardSubmissionScript(self):
        # If the remote computer is a cluster then you may need to submit jobs through a queuing system. This will act as a template for a standard job submission script but with no actual code to execute. The code to execute will be passed to the function as a list so that the child class doesn't have to keep being re-written everytime a new type of job is performed on the cluster.

    @abstractmethod
    def checkQueue(self):
        # Often a program needs to monitor the job queue of a cluster, queuing systems vary so this is left as an abstract method.
        pass

    @abstractmethod
    def checkDiskUsage(self):
        # To avoid running out of disk space on the remote computer one can call this function to find out how much is used and how much is available. Whilst this is fairly standard across *NIX home computers clusters and any kind of computer with multi-user storage often has custom functions for this and so is left as an abstract function.
        pass

    # INSTANCE METHODS
    def transferFile(self, source, destination, rsync_flags = "-aP"):
        """
        Uses rsync with specified flags (unspecified uses "-aP") to send a file to the remote computer. Source and destination only need the actual paths as the SSH connection will be done automatically.

        IMPORTANT: If you are unfamiliar with rsync then please read the manual or an introductory guide before using this function. Please remeber that a directory ending with no forward-slash indicates to copy the the directory and a directory ending with a forward-slash indicates to copy just the contents of the directory.

        Args:
            source (str): path and filename of the file to transfer.
            destination (str): path to the destination directory.
            rsync_flags (str): any flags that you would like when copying. This defaults to -aP (a is archive mode i.e. -rlptgoD (see manual for more information) and P is --progress which produces a progress bar.

        Returns:
            output_dict (dict): returns the 'return_code' from subprocess.call(rsync_cmd, shell=True).
        """
        rsync_cmd = "rsync " + rsync_flags + " " + source + " " + self.ssh_config_alias + ":" + destination
        output = subprocess.call(rsync_cmd, shell=True)
        output_dict = {}
        output_dict['return_code'] = output

        return output_dict

    def sendCommand(self, list_of_shell_commands):
        """
        'sendCommand' and 'getOutput' are the two main commands to communicate with a remote computer. This is easier for more complex commands.Note that the SSH connection is created automatically using the class variables and so does not need to be in the list_of_shell_commands.

        BEWARE: 'list_of_shell_commands' could be used as a vector for malicious injection and so in certain projects (specifically when untrusted users have direct access to this) could pose significant security risks.

        Args:
            list_of_shell_commands (list of strings): Each string in the list should be one whole command that will run in the shell language that runs on the remote computer. This list should NOT contain the SSH connection command. For example ['mkdir -p new_dir', 'cd new_dir'] will create the new_dir directory on the remote machine and the move into it.

        Returns:
            output_dict (dict): Has keys 'return_code', 'stdout', and 'stderr'.
        """

        # the -T flag in ssh is there because if you don't it opens a new instance of ssh everytime this function is run. It doesn't take long until your computer reaches it's maximum processes and then suddenly nothing can do anything because all the possible processes are being taken up by ssh instances not doing anything.
        sshProcess = subprocess.Popen(['ssh', '-T', self.ssh_config_alias], stdin=subprocess.PIPE, stdout = subprocess.PIPE, universal_newlines=True, bufsize=0)
        command = '\n'.join(list_of_shell_commands)
        print("command = ", command)
        out, err = sshProcess.communicate(command)
        return_code =  sshProcess.returncode
        sshProcess.stdin.close()
        output_dict = {}
        output_dict['return_code'] = return_code
        output_dict['stdout'] = out
        output_dict['stderr'] = err

        return output_dict

    # STATIC METHODS - I made these all static methods because I thought it might be handy to be able to use them without creating an instance.
    @staticmethod
    def checkSuccess(function, *args):
        """
        This function takes a function that requires a remote connection and makes sure that the actual command completes. If the connection can't be made then it it keeps trying whilst avoiding fast repeated connection attempts. After a connection can't be made for a whole day it trys once every 12 hours.
        
        Args:
            function (function): A function that makes a connection to a remote computer. This function MUST return a dictionary with atleast one element. This element have the key 'return_code' which returns the return code from the connection to the remote computer.
            *args (any combination of arguements): These will be the arguements needed to pass to function.
        Returns:
            output (unknown): Whatever function returns is saved as output and is returned.
            """

        # a list of the wait times (in seconds) between each loop should the connection keep failing
        # accumulative time:     15s 30s 45s  1m   6m  16m  30m   1hr   2hr   4hr    8hr   16hr   1day   2day   3day   4day   5day  6 day  7day
        wait_times = (3, 3, 3, 3, 3, 15, 15, 15, 300, 600, 840, 1800, 3600, 7200, 14400, 28800, 28800) #, 86400, 86400, 86400, 86400, 86400, 86400)

        # set flag to no successful connection attempt (successful exit code = 0)
        connection_success = 13
        for wait in wait_times:
            if connection_success != 0:
                try:
                    output = function(*args)
                    connection_success = output['return_code']
                except:
                    connection_success = 13

            else:
                break

            if connection_success !=0:
                print('Connection failed. Waiting ' + str(wait) + ' seconds before attempting to reconnect.')
                time.sleep(wait)

        # depending on the result either output the data or stop the simulation
        while connection_success != 0:
            43200
            try:
                output = function(*args)
                connection_success = output['return_code']
            except:
                connection_success = 13

            if connection_success !=0:
                print('Connection failed. Waiting ' + str(43200) + ' seconds before attempting to reconnect.')
                time.sleep(43200)

        return output

    @staticmethod
    def getOutput(commands_as_a_list):
        """
        This takes a list of commands (in the subprocess module style shell=False) and sends them via subprocess.check_output(commands_as_a_list). Note that if the commands_as_a_list does not start with an SSH command then the commands will be performed in a local shell.

        The shell=False is Python's prefered method to avoid malicious code injections.

        Args:
            commands_as_a_list (list of strings): Each contains a command, flag or option. For example ['mkdir', '-p', 'new_dir'] will create a new directory called new_dir on the local computer.
        Returns:
            raw_output (binary string): Returns the stdout as a binary string.
        """
        # if the return code is zero then subprocess returns the output of the command or otherwise raises an exception. We want to keep trying if there is an exception and so use the following code.
        try:
            raw_output = [0, subprocess.check_output(commands_as_a_list)]
        except Exception:
            raw_output = [1, None]
        return raw_output

