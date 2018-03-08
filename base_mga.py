from abc import ABCMeta, abstractmethod # Used to define abstract classes which is needed for the MGA class
import random # Some of the instance methods from the WholeCellModelBase class need the random library

class MGA(metaclass=ABCMeta):
    """
    This is an abstract class that all multi-generation algorithms inherit from. The idea is that the multi generation algorithm will want to repeatedly send jobs to (potentially multiple) computer clusters. With each generation the previous generation is made up of the parentjobs and the new generation is made up of the child jobs. This is so that information can be passed from generation to generation should the algorithm wish to do so.
    
    This class will assume that all connections are child classes of the base_connection.Connection class and all job submissions and job submission management classes are children of the relavent base_cluster_submissions class.
    """
    def __init__(self, dict_of_cluster_instances, MGA_name, relative2clusterBasePath_simulation_output_path, repetitions_of_a_unique_simulation):
        """
        This creates a basis for a multi-generation algorithm class.

        Args:
            dict_of_cluster_instances (dict): Each cluster instance is an object that is a child class of the base_connection.Connection class and represents a portal to a cluster or remote computer. The cluster instances are the values of the dictionary and the corresponding key is some string that represents some label of the cluster (these should be unique and it would be good for your records if you are consistent).
            MGA_name (str): The name of the multi-generation algorithm to be used as labels and names and records etc.
            relative2clusterBasePath_simulation_output_path (str): Each cluster connection instance has a base path depending on the cluster. This is the base path (i.e. the initial directory of this multi-generation algorithm) that will be apended to the cluster base path. More sub-directories will be created for each generation etc.
            repetitions_of_a_unique_simulation (int): The number of times each simulation needs to be repeated.
        """
        self.cluster_instances_dict = dict_of_cluster_instances
        self.generation_counter = None
        self.MGA_name = MGA_name
        self.relative2clusterBasePath_simulation_output_path = relative2clusterBasePath_simulation_output_path
        self.reps_of_unique_sim = repetitions_of_a_unique_simulation

    # instance methods
    def run(self):
        if self.generation_counter == None:
            self.generation_counter = 0
        while self.checkStop() != True:
            self.run_sim_out = self.runSimulations() 
            self.generation_counter += 1
            print('Next generation is ', self.generation_counter)

    def runSimulations(self):
        # get the new children
        child_name_to_set_dict = self.getNewGeneration()
        child_names_list = list(child_name_to_set_dict.keys())
        no_of_children = len(child_names_list)
        # split them across clusters
        no_of_clusters = len(self.cluster_instances_dict.keys())
        kos_per_cluster = int(no_of_children/no_of_clusters)
        # calculate the remainder
        remainder = no_of_children - (no_of_clusters * kos_per_cluster)
        child_name_to_set_dict_per_cluster = {}
        previous_idx = 0
        list_of_cluster_instance_keys = list(self.cluster_instances_dict.keys()) # we create this list because in most version of Python dictionary keys have no order so here we force a consistent order.
        for cluster_idx in range(len(self.cluster_instances_dict.keys())):
            cluster = list_of_cluster_instance_keys[cluster_idx]
            last_name_idx = (cluster_idx + 1) * kos_per_cluster
            if remainder > 0:
                last_name_idx += 1
                remainder -= 1

            child_name_to_set_dict_per_cluster[cluster] = {child_names_list[idx]: child_name_to_set_dict[child_names_list[idx]] for idx in range(previous_idx, last_name_idx)}
            previous_idx = last_name_idx
        
        # submit generation to the cluster
        dict_of_job_submission_insts = {}
        dict_of_job_management_insts = {}
        # create submission instances
        for cluster_connection in list_of_cluster_instance_keys:
            dict_of_job_submission_insts[cluster_connection] = self.createJobSubmissionInstance()

        # send all jobs to clusters 
        list_of_dict_of_job_management_instances = self.submitAndMonitorJobsOnCluster(dict_of_job_submission_insts)

        # convert list into the dict that the rest of the library is expecting
        dict_of_job_management_insts = {list_of_cluster_instance_keys[idx]: list_of_dict_of_job_management_instances[idx] for idx in range(len(list_of_dict_of_job_management_instances))}
        # update the fittest population
        for cluster_connection in list_of_cluster_instance_keys:
            self.updateFittestPopulation(dict_of_job_submission_insts[cluster_connection], dict_of_job_management_insts[cluster_connection])

        return

    # ABSTRACT METHODS
    @abstractmethod
    def submitAndMonitorJobsOnCluster(self, dict_of_job_submission_insts):
        # The job submission instance is an object that inherits from the base_cluster_submissions.BaseJobSubmission class. This takes a dictionary of job submission instance as an argument and then uses their methods to prepare and submit the jobs. From there this function can monitor the progress of the submission and do any other related work like process data and update databases etc. What needs to be done in this function needs to be defined at a higher level so this is left as an abstract method here.
        pass

    @abstractmethod
    def createJobSubmissionInstance(self):
        # The job submission instance is an object that inherits from the base_cluster_submissions.BaseJobSubmission class. This is left as an abstract method so that higher level programs can choose what kind of job submissions they want to create.
        pass

    @abstractmethod
    def getNewGeneration(self):
        pass

    @abstractmethod
    def getPopulationSize(self):
        pass

    @abstractmethod
    def mateTheFittest(self):
        pass

    @abstractmethod
    def getGenerationName(self):
        pass

    @abstractmethod
    def updateFittestPopulation(self):
        pass

    @abstractmethod
    def checkStop(self):
        pass

class WholeCellModelBase(MGA):
    def __init__(self, dict_of_cluster_instances, MGA_name, relative2clusterBasePath_simulation_output_path, reps_of_unique_sim):
        MGA.__init__(self, dict_of_cluster_instances, MGA_name, relative2clusterBasePath_simulation_output_path, relative2clusterBasePath_simulation_output_path)

    # All methods are static because they are common tools and it is useful to be able to access them wihtout creating an instance.

    # STATIC METHODS
    @staticmethod
    def random_combination(iterable_of_all_posible_indexs, length_of_combination):
            "Random selection from itertools.combinations(iterable_of_all_posible_indexs, length_of_combination)"
            pool = tuple(iterable_of_all_posible_indexs)
            n = len(pool)
            indices = sorted(random.sample(range(n), length_of_combination))

            return tuple(pool[i] for i in indices)

    @staticmethod
    def random_pick(list_of_options, probabilities):
        x = random.uniform(0, 1)
        cumulative_probability = 0.0
        for idx in range(len(list_of_options)):
            cumulative_probability += probabilities[idx]
            if x < cumulative_probability: 
                break

        return list_of_options[idx]

    @staticmethod
    def getGeneCodesToIdDict(conn, tuple_of_gene_codes):
        """Creates a dictionary who's keys are Joshua Rees 358 gene codes and values are the gene ID acording to our database."""
        gene_code_to_id_dict = conn.db_connection.convertGeneCodeToId(tuple_of_gene_codes)

        return gene_code_to_id_dict

    @staticmethod
    def invertDictionary(input_dict):
        """This function takes a dictionary and inverts it (assuming it's one to one)."""
        inverse_dict = {v: k for k, v in input_dict.items()}

        return inverse_dict

    @staticmethod
    def createIdxToIdDict(code_to_id_dict):
        list_of_ids = list(code_to_id_dict.values())
        # sort them into ascending order (just because the order of dicts aren't alwayys preserved and so provided we are using the same JR genes to start with we can compare the indexs provided they are ordered in ascending order) maybe not neccessary but avoiding hard to find bug later on
        list_of_ids.sort()
        idx_to_id_dict = {idx: list_of_ids[idx] for idx in range(len(list_of_ids))}

        return idx_to_id_dict

    @staticmethod
    def convertIdxToGeneId(gene_indexs_list, index_to_id_dict):
        """
        """
        # test input is of the right form
        if not (type(gene_indexs_list) is list and type(index_to_id_dict) is dict):
            raise TypeError('gene_indexs_list must be a list (even if only one value!) and index_to_id_dict must be a dictionary. Here type(gene_indexs_list)=', type(gene_indexs_list), ' and type(index_to_id_dict)=', type(index_to_id_dict))

        gene_id_list = [index_to_id_dict[idx] for idx in gene_indexs_list]

        return gene_id_list

    @staticmethod
    def convertGeneIdToCode(gene_id_list):
        """
        """
        # test input is of the right form
        if not (type(gene_id_list) is list):
            raise TypeError('gene_id_list must be a list (even if only one value!). Here type(gene_indexs_list)=', type(gene_indexs_list))

        gene_id_list = [index_to_id_dict[idx] for idx in gene_id_list]

    @staticmethod
    def getJr358Genes():
        """The function returns the 358 genes that Joshua Rees classified for potential KOs."""
        return ('MG_001', 'MG_003', 'MG_004', 'MG_005', 'MG_006', 'MG_007', 'MG_008', 'MG_009', 'MG_012', 'MG_013', 'MG_014', 'MG_015', 'MG_019', 'MG_020', 'MG_021', 'MG_022', 'MG_023', 'MG_026', 'MG_027', 'MG_029', 'MG_030', 'MG_031', 'MG_033', 'MG_034', 'MG_035', 'MG_036', 'MG_037', 'MG_038', 'MG_039', 'MG_040', 'MG_041', 'MG_042', 'MG_043', 'MG_044', 'MG_045', 'MG_046', 'MG_047', 'MG_048', 'MG_049', 'MG_050', 'MG_051', 'MG_052', 'MG_053', 'MG_055', 'MG_473', 'MG_058', 'MG_059', 'MG_061', 'MG_062', 'MG_063', 'MG_064', 'MG_065', 'MG_066', 'MG_069', 'MG_070', 'MG_071', 'MG_072', 'MG_073', 'MG_075', 'MG_077', 'MG_078', 'MG_079', 'MG_080', 'MG_081', 'MG_082', 'MG_083', 'MG_084', 'MG_085', 'MG_086', 'MG_087', 'MG_088', 'MG_089', 'MG_090', 'MG_091', 'MG_092', 'MG_093', 'MG_094', 'MG_097', 'MG_098', 'MG_099', 'MG_100', 'MG_101', 'MG_102', 'MG_476', 'MG_104', 'MG_105', 'MG_106', 'MG_107', 'MG_109', 'MG_110', 'MG_111', 'MG_112', 'MG_113', 'MG_114', 'MG_118', 'MG_119', 'MG_120', 'MG_121', 'MG_122', 'MG_123', 'MG_124', 'MG_126', 'MG_127', 'MG_128', 'MG_130', 'MG_132', 'MG_136', 'MG_137', 'MG_139', 'MG_141', 'MG_142', 'MG_143', 'MG_145', 'MG_149', 'MG_150', 'MG_151', 'MG_152', 'MG_153', 'MG_154', 'MG_155', 'MG_156', 'MG_157', 'MG_158', 'MG_159', 'MG_160', 'MG_161', 'MG_162', 'MG_163', 'MG_164', 'MG_165', 'MG_166', 'MG_167', 'MG_168', 'MG_169', 'MG_170', 'MG_171', 'MG_172', 'MG_173', 'MG_174', 'MG_175', 'MG_176', 'MG_177', 'MG_178', 'MG_179', 'MG_180', 'MG_181', 'MG_182', 'MG_183', 'MG_184', 'MG_186', 'MG_187', 'MG_188', 'MG_189', 'MG_190', 'MG_191', 'MG_192', 'MG_194', 'MG_195', 'MG_196', 'MG_197', 'MG_198', 'MG_200', 'MG_201', 'MG_203', 'MG_204', 'MG_205', 'MG_206', 'MG_208', 'MG_209', 'MG_210', 'MG_481', 'MG_482', 'MG_212', 'MG_213', 'MG_214', 'MG_215', 'MG_216', 'MG_217', 'MG_218', 'MG_221', 'MG_224', 'MG_225', 'MG_226', 'MG_227', 'MG_228', 'MG_229', 'MG_230', 'MG_231', 'MG_232', 'MG_234', 'MG_235', 'MG_236', 'MG_238', 'MG_239', 'MG_240', 'MG_244', 'MG_245', 'MG_249', 'MG_250', 'MG_251', 'MG_252', 'MG_253', 'MG_254', 'MG_257', 'MG_258', 'MG_259', 'MG_261', 'MG_262', 'MG_498', 'MG_264', 'MG_265', 'MG_266', 'MG_270', 'MG_271', 'MG_272', 'MG_273', 'MG_274', 'MG_275', 'MG_276', 'MG_277', 'MG_278', 'MG_282', 'MG_283', 'MG_287', 'MG_288', 'MG_289', 'MG_290', 'MG_291', 'MG_292', 'MG_293', 'MG_295', 'MG_297', 'MG_298', 'MG_299', 'MG_300', 'MG_301', 'MG_302', 'MG_303', 'MG_304', 'MG_305', 'MG_309', 'MG_310', 'MG_311', 'MG_312', 'MG_315', 'MG_316', 'MG_317', 'MG_318', 'MG_321', 'MG_322', 'MG_323', 'MG_324', 'MG_325', 'MG_327', 'MG_329', 'MG_330', 'MG_333', 'MG_334', 'MG_335', 'MG_517', 'MG_336', 'MG_339', 'MG_340', 'MG_341', 'MG_342', 'MG_344', 'MG_345', 'MG_346', 'MG_347', 'MG_349', 'MG_351', 'MG_352', 'MG_353', 'MG_355', 'MG_356', 'MG_357', 'MG_358', 'MG_359', 'MG_361', 'MG_362', 'MG_363', 'MG_522', 'MG_365', 'MG_367', 'MG_368', 'MG_369', 'MG_370', 'MG_372', 'MG_375', 'MG_376', 'MG_378', 'MG_379', 'MG_380', 'MG_382', 'MG_383', 'MG_384', 'MG_385', 'MG_386', 'MG_387', 'MG_390', 'MG_391', 'MG_392', 'MG_393', 'MG_394', 'MG_396', 'MG_398', 'MG_399', 'MG_400', 'MG_401', 'MG_402', 'MG_403', 'MG_404', 'MG_405', 'MG_407', 'MG_408', 'MG_409', 'MG_410', 'MG_411', 'MG_412', 'MG_417', 'MG_418', 'MG_419', 'MG_421', 'MG_424', 'MG_425', 'MG_426', 'MG_427', 'MG_428', 'MG_429', 'MG_430', 'MG_431', 'MG_433', 'MG_434', 'MG_435', 'MG_437', 'MG_438', 'MG_442', 'MG_444', 'MG_445', 'MG_446', 'MG_447', 'MG_448', 'MG_451', 'MG_453', 'MG_454', 'MG_455', 'MG_457', 'MG_458', 'MG_460', 'MG_462', 'MG_463', 'MG_464', 'MG_465', 'MG_466', 'MG_467', 'MG_468', 'MG_526', 'MG_470')
