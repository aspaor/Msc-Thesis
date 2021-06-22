
from random import random
from scipy.spatial import distance
import msprime
import os
os.environ["OMP_NUM_THREADS"] = "10"
import numpy as np
import math

import itertools
from itertools import chain
from sklearn.decomposition import PCA
from scipy import stats
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import matplotlib
import argparse
import dendropy
from dendropy.calculate import treecompare
from datetime import datetime
import allel
import allel.stats
from operator import itemgetter
from allel.stats.diversity import mean_pairwise_difference, mean_pairwise_difference_between
from allel.util import asarray_ndim,check_dim0_aligned, ensure_dim1_aligned
from sklearn.manifold import TSNE
from decimal import Decimal


class SupplementaryParameter(): pass
class ScoresParam: pass
class Scores: pass
class ScoresLists: pass
class TypeDistances : pass
class SimVariables(): pass
simvar = SimVariables()

startTime = datetime.now()
parameters_list = []

def main():
    
    de_list = []
    
    
    for i in DE:
        de = i.split(' ')
        de_list.append(de)
    de_list = [de_list]
    

    
    Samples = simvar.samples
    N = [simvar.N]
    R = [simvar.R]
    Length = [simvar.length]
    mutation_rate = [simvar.mutation_rate]
    recombination_rate = [simvar.recombination_rate]
    Distance,number_of_pops = make_distance(simvar.dist)
    Distance = [Distance] 
    
    migr_mat_list= make_migration_matrix(simvar.mig,number_of_pops)
    mig = migr_mat_list
    
    config = dict()
    config = {'de':de_list,'sample_list':Samples,'n':N,'r':R,'length':Length,'mutation_rate':mutation_rate,'recombination_rate':recombination_rate,'distance':Distance,'mig':mig}

    
    config, supplement_param = new_configuration_of_params(config,number_of_pops) #inputfile

    

    
    samples_values =  config['sample_list']
    supplement_param.population_configurations = []
    supplement_param.generation_time=25                               

    for each_Ne in config['n']: #Stop script with simulation of less than 3 populations
        if len(each_Ne)< 3:raise Exception('Populations less than 3,is not valid input') 
    
    for i,each_group in enumerate(config['n']):
        pop_config_old = []
        for j,each_N in enumerate(each_group):
            pop_config_old.append(msprime.PopulationConfiguration(initial_size=each_N,growth_rate=config['r'][i][j]))
        supplement_param.population_configurations.append(pop_config_old)
    
    simulation(supplement_param,config,samples_values,PCS = None) #simulation for 1000 reps
    
   
def check_if_is_a_list(x,parameter):
    ## Takes an item and checks if its type is a list and returns it. Otherwise it returns a list type item #
    
    if type(x) == list:
        if len(x)>1:
            if parameter not in parameters_list:
                parameters_list.append(parameter)

        return x
    else:
        return [x]

def createFolder(directory):
    ## Creates a folder if it does not exist already #
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

def write_parameters_in_file(parameter,sample,config,each_migr_matrix,each_mutation_rate,each_length,each_recombination_rate,each_demographic_event):
    folder = "./"+str(prefix)+'/'
    createFolder(folder)
    return folder

def checking_events():
    if event == 'bottleneck':
        if event not in parameters_list:
            parameters_list.append('bottleneck')
    elif event == 'recent_merge':
        if event not in parameters_list:
            parameters_list.append('recent_merge')
    elif event == 'oldest_merge':
        if event not in parameters_list:
            parameters_list.append('oldest_merge')
    elif event =='bottle_recent_merge':
        if event not in parameters_list:
            parameters_list.append('recent_merge')
            parameters_list.append('bottleneck')
    elif event =='bottle_oldest_merge':
        if event not in parameters_list:
            parameters_list.append('oldest_merge')
            parameters_list.append('bottleneck')
    return parameters_list


def writing_distances_in_file(folder,rep,info_di):
    
    with open(folder+'/distances'+'.txt','a') as f:
        f.write('rep {}\n'.format(rep))
        f.write("fst scores between pops {}\n".format(info_di['fst']))
        f.write("\neuclideian distance of centers of PCA data{}\n".format(info_di['pca_eucl']))
        f.write("\n\neuclideian distance of centers of TSNE data{}\n\n".format(info_di['tsne_eucl']))
                                              
    f.close()
    

                                            
def simulation(supplement_param,config,samples_values,PCS = None):
#def simulation(generation_time,mutation_rate,length,recombination_rate,population_configurations,migration_matrix,demographic_events,prefix,N,r,Distance_before_pca,recent_merge_list,oldest_merge_list,pop_bottled_list,pop_bottled,bottle_st,peer_st,bottle_end_list,peer_end,merged_pops,samples_values,PCS):
    ## Takes parameters as input to make an msprime simulation 
    
    parameter = 0
    number_of_pop = len(config['n'][0])
    keys = ['pca_eucl','fst','tsne_eucl']
    scor_di = {k: [] for k in keys}
    
    simulation_parameter_combinations = itertools.product(
        check_if_is_a_list(samples_values,'sample'),
        check_if_is_a_list(config['mig'],'migration'),
        check_if_is_a_list(supplement_param.population_configurations,'pop_config'),
        check_if_is_a_list(config['mutation_rate'],'mut_rate'),
        check_if_is_a_list(config['length'],'length'),
        check_if_is_a_list(config['recombination_rate'],'recomb_rate'),
        check_if_is_a_list(supplement_param.demographic_events,'demo_event'),
        check_if_is_a_list(config['distance'],'distance')
        )      

   

    for sim_param_comb in simulation_parameter_combinations:
        sample = sim_param_comb[0]
        
        each_migr_matrix = sim_param_comb[1]
        each_population_configuration = sim_param_comb[2]
        each_mutation_rate = sim_param_comb[3]
        each_length = sim_param_comb[4]
        each_recombination_rate = sim_param_comb[5]
        each_demographic_event = sim_param_comb[6]
        each_old_distance = sim_param_comb[7]
        
        inside_startTime = datetime.now()
        checking_events()       
        parameter+=1
        dicts_centers = {}
        tsne_dicts_centers = {}
        
        if number_of_pop ==len(each_population_configuration) == len(each_migr_matrix) == len(each_old_distance):
            
            
            PC_list = []
            samples = []
            number_of_total_diploids = 0
            lst_of_pop_samples = []

            folder = write_parameters_in_file(parameter,sample,config,each_migr_matrix,each_mutation_rate,each_length,each_recombination_rate,each_demographic_event)
            
            for pop_index in range(number_of_pop):
                
                samples += [msprime.Sample(pop_index,0)]*sample 
                number_of_total_diploids += int(sample/2)  
                lst_of_pop_samples += [int(sample/2)] 
            
            SIMULATION= msprime.simulate(samples=samples, population_configurations=each_population_configuration, migration_matrix=each_migr_matrix, mutation_rate=each_mutation_rate ,length=each_length, recombination_rate=each_recombination_rate,demographic_events=each_demographic_event,num_replicates=100)
            fst_after_reps = []
            
            
            scores_lists = {l: [] for l in keys}
            
            for rep,each_simulation in enumerate(SIMULATION):
                
                print('\n rep ',rep,'\n-------')
                
                samples = sample #sample[0]
                genot_matrix = each_simulation.genotype_matrix()
                snps, number_of_samples = genot_matrix.shape
                print('SHAPE_DATA_BEFORE_FST ta snps:{} kai ta nsamples:{}'.format(snps,number_of_samples))
                
                genot_matrix = ld_clumping(genot_matrix)
                type_dist = TypeDistances()
                fst_after_reps,type_dist.fst_bet_pops = new_fst(genot_matrix,fst_after_reps,number_of_pop)
                
                
                X_transformed,PC_list,PC=PCA_plot(genot_matrix,number_of_pop,lst_of_pop_samples,PC_list) #[0][1]
                
                PC = choose_pc(PC,PCS)
                shrinked_matrix = keep_pcs(X_transformed,PC) #o pinakas meta ti pca pou kratisame osa PCs (columns) thelame
        
                X_tsne = TSNE_process(shrinked_matrix,lst_of_pop_samples)
                
                dicts_centers,type_dist.euclideian_dist_of_pops= find_centers(shrinked_matrix,number_of_pop,dicts_centers,each_old_distance) #P1,P2,P3 oi prwtes metablites
                tsne_dicts_centers,type_dist.tsne_euclideian_dist_of_pops = find_centers(X_tsne,number_of_pop,tsne_dicts_centers,each_old_distance) #P1,P2,P3 oi prwtes metablites

                info_di = {
                "pca_eucl": type_dist.euclideian_dist_of_pops,
                'fst':type_dist.fst_bet_pops,
                'tsne_eucl':type_dist.tsne_euclideian_dist_of_pops,
                }
                
                print(shrinked_matrix.shape,'shape auto pou emeine apo to pinaka\n')
                
                scores_lists = calculate_similarity(each_old_distance,info_di,number_of_pop,scores_lists)
                
                writing_distances_in_file(folder,rep,info_di)

             
            
            scor_di = computation_of_scores(scores_lists,scor_di)
            
            
            os.system("command -v osascript > /dev/null && osascript -e 'display notification \"Plot Done\"'") #sto allagi param htan prin
        #print('Time of each variable alteration ',datetime.now() - inside_startTime) 


    parameter_to_plot = plot_parameters_scores(parameters_list,samples_values,supplement_param,config)           
    #print(parameter_to_plot,'focused parameters')
    
    
    print('___SCORES___\n\n','tsne_euc',scor_di['tsne_eucl'])
    print('\npca_eucl',scor_di['pca_eucl'])
    print('\nfstscore',scor_di['fst'])

    
    if len(parameters_list) ==1:
        if parameter == 'bottleneck' or parameter == 'recent_merge' or parameter == 'oldest_merge':
            parameter_to_plot = [parameter_to_plot]
        
        if supplement_param.merged_pops != []:
            supplement_param.merged_pops = supplement_param.merged_pops[0]
        writing_scores(supplement_param,scor_di,parameters_list[0],parameter_to_plot,folder,samples_values) #writes fst,euclideian and nucl_diff scores in output files

    if len(parameters_list) == 2:
        if supplement_param.merged_pops != []:
            supplement_param.merged_pops = supplement_param.merged_pops[0]
            
            
        else:
            if 'merge' in event:
                supplement_param.merged_pops = supplement_param.merged_pops[0]
            
            


        writing_scores(supplement_param,scor_di,parameters_list,parameter_to_plot,folder,samples_values) #writes fst,euclideian and nucl_diff scores in output files

def choose_pc(PC,PCS):
    if PC<2: #gia times 0,1
        PC = 2
    PC = int(PC)
    if PCS != None:
        PC = PCS
    return PC

def calculate_similarity(each_old_distance,info_di,number_of_pop,scores_lists):
    ''' Calls similarity, fst_similarity, multiple_pops_similarity functions and fills scores_lists
    
        Parameters
        ----------
        each_old_distance: arraylike, contains distance values of populations before any process ,(scenario distance)
        info_di: dict,contains lists of different type of distances of populations(pca_eucl,fst,tsne_eucl etc)
        number_of_pop: number of populations
        scores_lists: dict,contains lists of number of success (0,1) from recognition process of populations
    '''
    
    
    for type_metric,metric_val in info_di.items():
        if number_of_pop == 3:
            if type_metric != 'fst':
                result_of_similarity = similarity(each_old_distance,metric_val,number_of_pop)
            else:
                result_of_similarity = fst_similarity(metric_val,each_old_distance,number_of_pop)

        elif number_of_pop > 3:
            result_of_similarity = multiple_pops_similarity(each_old_distance,metric_val,number_of_pop)

        scores_lists[type_metric].append(result_of_similarity)
    return scores_lists



def new_configuration_of_params(config,number_of_pops):
    
    ## Configuration of parameters using as input a txt file
    
    supplement_param = SupplementaryParameter()
    merge_time = 0
    supplement_param.recent_merge_list = []
    supplement_param.oldest_merge_list = []
    N_list = []
    r_list = []
    mut_rat_list = []
    length_list = []
    recombin_list = []
    samp_list = []
    distance_time_list = []
    migr_list = []
    event_list = []
    supplement_param.demographic_events = []
    distance_matrix_list = []
    
   
    for each_key in config.keys():
        if 'n' == each_key:
            for each_N,_ in enumerate(config['n']):
                N = [int(i) for i in config['n'][each_N]]
                N_list.append(N)
            N_list = sorted(N_list, key=itemgetter(0))
            
            config['n']= N_list
            
            number_pop = len(N)  
            
        elif 'r' == each_key:
            for each_r,_ in enumerate(config['r']):
                r = [float(i) for i in config['r'][each_r]]
                r_list.append(r)

            r_list = sorted(r_list, key=itemgetter(0))
            config['r']= r_list
            
        elif 'mutation_rate' == each_key:
            for each_mutation_rate,_ in enumerate(config['mutation_rate']):
                for i in config['mutation_rate'][each_mutation_rate]:
                    mut_rat_list.append(float(i))
                    mut_rat_list = sorted(mut_rat_list)
                    
            config['mutation_rate']= mut_rat_list
                
        elif 'length' == each_key:
            for each_length,_ in enumerate(config['length']):
                for i in config['length'][each_length]:
                    length_list.append(float(i))
                    length_list = sorted(length_list)
                    
            config['length']= length_list
            
        elif 'recombination_rate' == each_key:
            for each_recombination_rate,_ in enumerate(config['recombination_rate']):
                for i in config['recombination_rate'][each_recombination_rate]:
                    recombin_list.append(float(i))
                    recombin_list = sorted(recombin_list)
                    
            config['recombination_rate']= recombin_list
                
        elif 'sample_list' == each_key:
            
            for each_sample_list in config['sample_list']:
                samp_list +=[int(each_sample_list)]
            
            
            samp_list = sorted(samp_list)
            config['sample_list']= samp_list
            
            
        elif 'distance' == each_key:
            for each_distance,_ in enumerate(config['distance']):
                
                distance_time_list.append(config['distance'][each_distance])
            config['distance']= distance_time_list
            
            for i in config['distance']:
                distance_matrix_list.append( np.asarray(i,dtype=int).reshape(number_pop,number_pop)) #list with migration matrices
            config['distance']= distance_matrix_list
            
            
        elif 'mig' == each_key:
            
            for each_migration ,_ in enumerate(config['mig']):
                migration_matrix = [float(j) for i in config['mig'][each_migration] for j in i]
                migration_matrix = np.asarray(migration_matrix)
                migration_matrix = migration_matrix.tolist()
                migration_matrix = [migration_matrix]
                migr_list.append(migration_matrix)
                
            config['mig'] = migr_list
            migr_list = []
            for i in config['mig']:
                migration = list(chain(*i))
                migr_list.append( np.asarray(migration,dtype=float).reshape(number_pop,number_pop).tolist()) #list with migration matrices
            config['mig'] = migr_list
            
            
        elif 'de' == each_key:
            
            accepted_pop_indeces = range(0,number_of_pops)
            EVENT = []
            supplement_param.merged_pops = []
            supplement_param.pop_bottled_list = []
            count = 0
            supplement_param.pop_bottled = -1
            is_there_ppc = 0
            supplement_param.bottle_end_list = []

            for every_event_group in config[each_key]:
                event_list = []
                for every_event in every_event_group:
                    
                    
                    if every_event[1].lower() == 'ppc':
                        is_there_ppc +=1
                        if int(every_event[4]) in accepted_pop_indeces:
                            
                            if count==0: #first time that a pop changes
                                supplement_param.pop_bottled = int(every_event[4])
                                supplement_param.bottle_start = int(every_event[2])
                                
                                supplement_param.people_start = int(every_event[3])
                                count+=1

                            elif count!=0 and supplement_param.pop_bottled!=-1:
                                if supplement_param.pop_bottled==int(every_event[4]):
                                    bottle_end = int(every_event[2])
                                    supplement_param.pop_bottled_list.append(supplement_param.pop_bottled)
                                    count = 0
                                    supplement_param.pop_bottled =-1
                                    if supplement_param.bottle_start != bottle_end:
                                        supplement_param.bottle_end_list.append(bottle_end)
                                        supplement_param.people_end = int(every_event[3])
                                    
                            event_list.append((every_event[1],int(every_event[2]),int(every_event[3]),int(every_event[4])))
                        else:
                            raise Exception('wrong index of population')

                    if is_there_ppc == 0: #if no ppc there will no be any bottleneck
                        
                        supplement_param.bottle_start = 0
                        supplement_param.bottle_end_list = []
                        supplement_param.people_start = 0
                        supplement_param.people_end = 0

                        if every_event[1].lower() == 'mm':
                            if int(every_event[4]) in accepted_pop_indeces and int(every_event[3]) in accepted_pop_indeces:
                                event_list.append((every_event[1],int(every_event[2]),int(every_event[3]),int(every_event[4]),float(every_event[5])))
                            else:
                                raise Exception('wrooong index of population')
                        elif every_event[1].lower() == 'mrc':
                            if int(every_event[4]) in accepted_pop_indeces and int(every_event[5]) in accepted_pop_indeces :
                                event_list.append((every_event[1],int(every_event[2]),float(every_event[3]),int(every_event[4]),int(every_event[5])))
                            
                            else:
                                raise Exception('Wrong index of population')
                    
                    if every_event[1].lower() == 'mm':
                        if event == "recent_merge" or event =='bottle_recent_merge':
                            if merge_time == 0 :
                                merge_time = int(every_event[2])
                                supplement_param.merged_pops.append(((int(every_event[3])+1),int(every_event[4])+1))
                            else:
                                if merge_time > int(every_event[2]):
                                    merge_time = int(every_event[2])
                        if event == "oldest_merge" or event == 'bottle_oldest_merge':
                            if merge_time == 0 :
                                merge_time = int(every_event[2])
                            else:
                                if merge_time < int(every_event[2]):
                                    merge_time = int(every_event[2])
                                    supplement_param.merged_pops.append(((int(every_event[3])+1),int(every_event[4])+1))
                            
                        if int(every_event[4]) in accepted_pop_indeces and int(every_event[3]) in accepted_pop_indeces:
                            event_list.append((every_event[1],int(every_event[2]),int(every_event[3]),int(every_event[4]),float(every_event[5])))
                        else:
                            raise Exception('wrooong index of population')
                
                    if every_event[1].lower() == 'mrc':
                        if int(every_event[4]) in accepted_pop_indeces and int(every_event[5]) in accepted_pop_indeces :
                            event_list.append((every_event[1],int(every_event[2]),float(every_event[3]),int(every_event[4]),int(every_event[5])))
                        else:
                            raise Exception('Wrong index of population')
                        
                sort_events = sorted(event_list, key=lambda x: x[1])
                EVENT.append(sort_events)
                if (event == 'recent_merge') or (event == 'bottle_recent_merge'):
                    supplement_param.recent_merge_list.append(merge_time)
                    
                    
                elif (event == 'oldest_merge') or (event == 'bottle_oldest_merge'):
                    supplement_param.oldest_merge_list.append(merge_time)
            
            
            demo_ev = []
            for every in EVENT:
                demo_ev =[]
                for each in every:
                    if each[0] == 'ppc':
                        demo_ev.append( msprime.PopulationParametersChange(time=each[1]/25,initial_size=each[2],population_id=each[3]))
                    elif each[0] == 'mm':
                        demo_ev.append( msprime.MassMigration(time=each[1]/25,source=each[2],destination=each[3],proportion = each[4]))

                    elif each[0] == 'mrc':
                        demo_ev.append( msprime.MigrationRateChange(time=each[1]/25 , rate=each[2], matrix_index=(each[3],each[4])))

                supplement_param.demographic_events.append(demo_ev)
            
        elif each_key == ' ' or each_key == '\n' or each_key == '' or each_key =='  ':
            next
        else:
            print(each_key,'got in this Error case')
            raise Exception('not valid parameter')
    #print('---------------ta configssss mouuuuuu')
    # for i in config.items():
    #     print(i,'wwwwwcwcwcwc')
    #print('------f--------\n')
    
    return config,supplement_param





def configuration_of_params(file):
    #Configuration of parameters using as input a txt file
    
    with open(file,'r') as f: 
        li = 1
        config = {}
        bottleneck_commands = []
        event_creation = [] 
        all_events = []
        
        multippc=False
        for line in f:
            if line == "\n": 
                li+=1
                
                continue
            else: #if has characters
                # cmd_arguments = ['de', 'ppc', '2000', '100', '2']
                cmd_arguments = line.lower().strip().rsplit(' ')

                name = cmd_arguments[0].lower()
                val = cmd_arguments[1:]
                
                if (name in config.keys()) and (name !='de'):
                    config[name].append(val)
                elif (name not in config.keys()) and (name !='de'):
                    config[name]=[val]
                elif name == 'de':
                    if '-' not in cmd_arguments[2]:
                        event_creation.append(cmd_arguments)
                        all_events.append(cmd_arguments)
                    else:
                        multippc=True

                        ppc_times = []
                        prin,meta = cmd_arguments[2].split('-')
                        
                        for t in range(int(prin),int(meta),10000):
                            ppc_times.append(t)
                            
                        # Replace the time in current command, then append it
                        for time in ppc_times:
                            
                            tmp = cmd_arguments[:]
                            tmp[2] = time

                            bottleneck_commands.append(tmp)
            
        if multippc:
            all_events = []
            for bcommand in bottleneck_commands: # ppc1 10.000 10.700 
                events = []
                for x in event_creation:
                    events.append(x)
                events.append(bcommand)
                all_events.append(events)
            
            config[name]=all_events
        else:
            config[name]=[all_events]
    

    supplement_param = SupplementaryParameter()
    merge_time = 0
    supplement_param.recent_merge_list = []
    supplement_param.oldest_merge_list = []
    N_list = []
    r_list = []
    mut_rat_list = []
    length_list = []
    recombin_list = []
    samp_list = []
    distance_time_list = []
    migr_list = []
    event_list = []
    supplement_param.demographic_events = []
    distance_matrix_list = []

    
    
    for each_key in config.keys():
        if 'n' == each_key:
            for each_N,_ in enumerate(config['n']):
                N = [int(i) for i in config['n'][each_N]]
                N_list.append(N)
            N_list = sorted(N_list, key=itemgetter(0))
            
            config['n']= N_list
            #print(config['n'],'FTIAKSIK')
            
            number_pop = len(N)  
            
        elif 'r' == each_key:
            for each_r,_ in enumerate(config['r']):
                r = [float(i) for i in config['r'][each_r]]
                r_list.append(r)

            r_list = sorted(r_list, key=itemgetter(0))
            config['r']= r_list
            
        elif 'mutation_rate' == each_key:
            for each_mutation_rate,_ in enumerate(config['mutation_rate']):
                for i in config['mutation_rate'][each_mutation_rate]:
                    mut_rat_list.append(float(i))
                    mut_rat_list = sorted(mut_rat_list)
                    
            config['mutation_rate']= mut_rat_list
                
        elif 'length' == each_key:
            for each_length,_ in enumerate(config['length']):
                for i in config['length'][each_length]:
                    length_list.append(float(i))
                    length_list = sorted(length_list)
                    
            config['length']= length_list
            
        elif 'recombination_rate' == each_key:
            for each_recombination_rate,_ in enumerate(config['recombination_rate']):
                for i in config['recombination_rate'][each_recombination_rate]:
                    recombin_list.append(float(i))
                    recombin_list = sorted(recombin_list)
                    
            config['recombination_rate']= recombin_list
                
        elif 'sample_list' == each_key:
            for each_sample_list in config['sample_list']:
                sample_list = [int(each_sample_list[0])]
                samp_list.append(sample_list)
            samp_list = sorted(samp_list, key=itemgetter(0))
            config['sample_list']= samp_list
               
        elif 'distance' == each_key:
            for each_distance,_ in enumerate(config['distance']):
                Distance = [int(i) for i in config['distance'][each_distance]]
                distance_time_list.append(Distance)
            config['distance']= distance_time_list
            
            for i in config['distance']:
                distance_matrix_list.append( np.asarray(i,dtype=int).reshape(number_pop,number_pop)) #list with migration matrices
            config['distance']= distance_matrix_list

        elif 'mig' == each_key:
            
            for each_migration ,_ in enumerate(config['mig']):
                migration_matrix = [float(i) for i in config['mig'][each_migration]]
                migration_matrix = np.asarray(migration_matrix)
                migration_matrix = migration_matrix.tolist()
                migration_matrix = [migration_matrix]
                migr_list.append(migration_matrix)
            config['mig'] = migr_list
            migr_list = []
            for i in config['mig']:
                migration = list(chain(*i))
                migr_list.append( np.asarray(migration,dtype=float).reshape(number_pop,number_pop).tolist()) #list with migration matrices
            config['mig'] = migr_list
            
        elif 'DE'.lower() == each_key:
            
            accepted_pop_indeces = range(0,number_pop)
            EVENT = []
            supplement_param.merged_pops = []
            supplement_param.pop_bottled_list = []
            count = 0
            supplement_param.pop_bottled = -1
            is_there_ppc = 0
            supplement_param.bottle_end_list = []

            for every_event_group in config[each_key]:
                event_list = []
                
                for every_event in every_event_group:
                    
                    
                    if every_event[1].lower() == 'ppc':
                        is_there_ppc +=1
                        if int(every_event[4]) in accepted_pop_indeces:
                            
                            if count==0: #first time that a pop changes
                                supplement_param.pop_bottled = int(every_event[4])
                                supplement_param.bottle_start = int(every_event[2])
                                
                                supplement_param.people_start = int(every_event[3])
                                count+=1

                            elif count!=0 and supplement_param.pop_bottled!=-1:
                                if supplement_param.pop_bottled==int(every_event[4]):
                                    bottle_end = int(every_event[2])
                                    supplement_param.pop_bottled_list.append(supplement_param.pop_bottled)
                                    count = 0
                                    supplement_param.pop_bottled =-1
                                    if supplement_param.bottle_start != bottle_end:
                                        supplement_param.bottle_end_list.append(bottle_end)
                                        supplement_param.people_end = int(every_event[3])
                                    
                            event_list.append((every_event[1],int(every_event[2]),int(every_event[3]),int(every_event[4])))
                        else:
                            raise Exception('wrong index of population')

                    if is_there_ppc == 0: #if no ppc there will no be any bottleneck
                        
                        supplement_param.bottle_start = 0
                        supplement_param.bottle_end_list = []
                        supplement_param.people_start = 0
                        supplement_param.people_end = 0

                        if every_event[1].lower() == 'mm':
                            if int(every_event[4]) in accepted_pop_indeces and int(every_event[3]) in accepted_pop_indeces:
                                event_list.append((every_event[1],int(every_event[2]),int(every_event[3]),int(every_event[4]),float(every_event[5])))
                            else:
                                raise Exception('wrooong index of population')
                        elif every_event[1].lower() == 'mrc':
                            if int(every_event[4]) in accepted_pop_indeces and int(every_event[5]) in accepted_pop_indeces :
                                event_list.append((every_event[1],int(every_event[2]),float(every_event[3]),int(every_event[4]),int(every_event[5])))
                            
                            else:
                                raise Exception('Wrong index of population')
                    
                    if every_event[1].lower() == 'mm':
                        if event == "recent_merge" or event =='bottle_recent_merge':
                            if merge_time == 0 :
                                merge_time = int(every_event[2])
                                supplement_param.merged_pops.append(((int(every_event[3])+1),int(every_event[4])+1))
                            else:
                                if merge_time > int(every_event[2]):
                                    merge_time = int(every_event[2])
                        if event == "oldest_merge" or event == 'bottle_oldest_merge':
                            if merge_time == 0 :
                                merge_time = int(every_event[2])
                            else:
                                if merge_time < int(every_event[2]):
                                    merge_time = int(every_event[2])
                                    supplement_param.merged_pops.append(((int(every_event[3])+1),int(every_event[4])+1))
                            
                        if int(every_event[4]) in accepted_pop_indeces and int(every_event[3]) in accepted_pop_indeces:
                            event_list.append((every_event[1],int(every_event[2]),int(every_event[3]),int(every_event[4]),float(every_event[5])))
                        else:
                            raise Exception('wrooong index of population')
                
                    if every_event[1].lower() == 'mrc':
                        if int(every_event[4]) in accepted_pop_indeces and int(every_event[5]) in accepted_pop_indeces :
                            event_list.append((every_event[1],int(every_event[2]),float(every_event[3]),int(every_event[4]),int(every_event[5])))
                        else:
                            raise Exception('Wrong index of population')
                        
                sort_events = sorted(event_list, key=lambda x: x[1])
                EVENT.append(sort_events)
                if (event == 'recent_merge') or (event == 'bottle_recent_merge'):
                    supplement_param.recent_merge_list.append(merge_time)
                    
                    
                elif (event == 'oldest_merge') or (event == 'bottle_oldest_merge'):
                    supplement_param.oldest_merge_list.append(merge_time)
            
            
            demo_ev = []
            for every in EVENT:
                demo_ev =[]
                for each in every:
                    if each[0] == 'ppc':
                        demo_ev.append( msprime.PopulationParametersChange(time=each[1]/25,initial_size=each[2],population_id=each[3]))
                    elif each[0] == 'mm':
                        demo_ev.append( msprime.MassMigration(time=each[1]/25,source=each[2],destination=each[3],proportion = each[4]))

                    elif each[0] == 'mrc':
                        demo_ev.append( msprime.MigrationRateChange(time=each[1]/25 , rate=each[2], matrix_index=(each[3],each[4])))

                supplement_param.demographic_events.append(demo_ev)
            
        elif each_key == ' ' or each_key == '\n' or each_key == '':
            next
        else:
            print(each_key,"\n This error came up")
            raise Exception('not valid parameter')

    #print('---------------ta configssss mouuuuuu')
    # for i in config.items():
    #     print(i)
    print('-----------------------------------------------------------------------------------------------------------------\n')
    
    return config,supplement_param




def TSNE_process(ld_clumped_matrix,lst_of_pop_samples):

    diploid_per_pop = lst_of_pop_samples[0]
    
    X_tsne_reduced=TSNE(perplexity=40).fit_transform(ld_clumped_matrix)
    print('Xtsne',X_tsne_reduced.shape)
    
    return X_tsne_reduced
    
def PCA_plot(ld_clumped_matrix,number_of_populations,lst_of_pop_samples,PC_list):
    ## Implements PCA to the simulated data in order to reduce dimensionality. Returns pca-transformed data, list of PCS and ideal number of PCs for each run
    
    # X = np.array(converted_list).reshape((count,diploids))   #DxN   features xsamples etsi ta dinei to vcf 
    # print('arxiko shape of diploid_data before pca(snps,diploid) ',X.shape)
    
    ld_clumped_matrix = ld_clumped_matrix.T
    # create the labels for pca plot
    population_li = [i+1 for i in range(number_of_populations)]
    label = [np.repeat(population_li[i],lst_of_pop_samples[i]) for i in range(len(population_li))] #labels repeat same label for times = sample for each pop
    Y=np.asarray(list(itertools.chain.from_iterable(label)))
    
    scaler = StandardScaler()
    X_std=scaler.fit_transform(ld_clumped_matrix) #standarization of values
    pca=PCA()
    X_std=pca.fit_transform(X_std)
    X_transformed = X_std
    explained_variance= np.round(pca.explained_variance_ratio_ *100,decimals = 1)
    
    idx1 = np.where(np.abs(np.diff(explained_variance))<=0.000001)[0] # no difference in explained variance, points where slope of differentiation doesn't reach zero,doesn't change so much
    
    if len(idx1)!=0:
        idx2 = np.where(explained_variance<=0.15*np.max(explained_variance))[0] # where explained is below a threshold, % of PC1 explained var, where explained var is 4% of maximum explained var
        
        if len(idx2)!=0:
            if len(set(idx1).intersection(set(idx2)))==0: # if idx1, idx2 don't have common numbers choose first elemnt of idx1,idx1[0] chosen
                
                PC = idx1[0]
            else:
                # Find the common idx between those two, and subtract 1
                idx = np.intersect1d(idx1, idx2)[0]
                PC = idx-1 # prinicipal components
            
        else: 
            PC=len(explained_variance) #otherwise take all prinicipal components
    else: 
        PC=len(explained_variance) # otherwise prinicipal components


    #plt.figure() ; plt.plot(explained_variance); plt.plot(explained_variance, '.'); plt.show()#plt.plot(PC, explained_variance[PC], '*') ;plt.show()
    
    
    
    
    PC_list.append(PC)
    

    
    return X_transformed,PC_list,PC



def keep_pcs(matrix,PCS): 
    ## Takes as input the number of PCs wanted to be kept and the matrix after pca. Returns a matrix with columns equal to the wanted number of PCs 

    _,columns = matrix.shape
    if PCS <= columns:
        mat = matrix[:,0:PCS]
        return mat
    else:
        raise Exception("PCs are more than matrix's columns,try again")



def find_centers(pin,numb_of_pop,dicts_centers,each_old_distance): 
    ## Finds centers of each pop, by calculating mean x mean y
    
    pops= np.split(pin,numb_of_pop) #kaneis split ton x se 3 pinakes (populations)
    group = 1
    list_of_centers = []
    dictionary = {}
    
    for each_arr in pops: #gia kathe plithismo
        center_of_group = list(np.mean(each_arr, axis=0)) #find mean x mean y,so center of each pop is (mean_x,mean_y)
        list_of_centers.append(center_of_group) #to bazeis sti lista me ta kentra
        group+=1
    
    euclideian_dist_of_pops = distance_bet_pop(list_of_centers)
    
    if dicts_centers == {}:
        keys = range(group-1) #kleidia tou dict oi plthusmoi 0,1,2 
        for i in keys:
            dicts_centers[i] = list_of_centers[i] #values of old pops' centers ,in the end it will contain all centers after reps
    else:
        keys = range(group-1)
        for i in keys:
            dictionary[i] = list_of_centers[i] #values of pops' new centers
        for key,val in dictionary.items():
            if key in dicts_centers:
                dicts_centers[key] += val 
        
    return dicts_centers,euclideian_dist_of_pops


def finding_mean_centers(dicts_centers,reps,PCS,number_of_pop): 
    ## Takes a dictionary with keys number of pop and values a list of centers of each pop after reps and calculates the mean center of each population
    
    plithsumo = 1
    mean_centers_list = []
    centers_of_pops_list = []
    for center in dicts_centers.values():
        center = np.asarray(center).reshape(reps,PCS) #center_group with centers of each population
        centers_of_pops_list.append(center)  
        center_of_pop = tuple(np.mean(center,axis=0))
        mean_centers_list.append(center_of_pop) #mean_center_list involves mean center of each pop
        plithsumo+=1
    
    centers= np.split(np.asarray(centers_of_pops_list).reshape(number_of_pop*reps,PCS),number_of_pop) #C1,C2,C3 are mean centers of pops
    
    return mean_centers_list,centers

def distance_bet_pop(lista):
    ## Takes mean_center_list and finds euklideian distance between mean_centers of pops d(1,2),d(1,3),d(2,3)
    euclideian_dist_of_pops = distance.pdist(lista) #pairwise calculation
    
    return euclideian_dist_of_pops


def similarity(each_old_distance,distances_of_pops,number_of_pop):
    ## Compares 2 matrices(distances of pops before and after PCA)
    
    sumofel_before = 0
    sumofel_after = 0
    elements1 = 0
    elements2 = 0

    anw_before = each_old_distance[np.triu_indices(3, k = 1)]
    anw_after = distances_of_pops
    for i in anw_before :
        if i!=min(anw_before):
            sumofel_before = sumofel_before +i
        if i==min(anw_before):
            elements1+=1

    if elements1 == 2:
        print('there are elements with same value')
        return 0 
    else:
        mean_of_anw_before = sumofel_before/2 #find mean of 2 biggest elements
     
    for i in anw_after :
        if i!=min(anw_after):
            sumofel_after = sumofel_after +i
        if i==min(anw_after):
            elements2+=1

    if elements2 == 2 or elements2==3: ##CHANGED
        print('there are elements with same value')
        return 0 
    else:
        mean_of_anw_after = sumofel_after/2 #find mean of 2 biggest elements

    minof_before = min(anw_before)
    trimmed_list_before = [x for x in anw_before if x!=minof_before]  
    diafora_trimmed_before = abs(trimmed_list_before[0]-trimmed_list_before[1])
    diafora_small_numbers_before = abs(min(trimmed_list_before)-minof_before)

    if diafora_small_numbers_before < diafora_trimmed_before:
        return 0 

    minof_after = min(anw_after)
    trimmed_list_after = [x for x in anw_after if x!=minof_after]  
    diafora_trimmed_after = abs(trimmed_list_after[0]-trimmed_list_after[1])
    diafora_small_numbers_after = abs(min(trimmed_list_after)-minof_after)
    
    if diafora_small_numbers_after < diafora_trimmed_after:
        return 0

    if (mean_of_anw_before > min(anw_before)) and  (mean_of_anw_after > min(anw_after)) :
        if np.argmin(anw_before) == np.argmin(anw_after):
            return 1
        else:
            return 0 
    else:
        return 0   


def create_tree(lista,npop):
    ## Takes a list as input with numbers in time and creates a tree

    string =''
    i = 1
    pop = 1
    thesi = 1
    for i in range(1,len(lista)+1):
        if i == 1:
            new_string = str(pop)+':'+str(sum(lista[0:npop-pop])/len(lista[0:npop-pop]))
            thesi = thesi + npop-2
            pop+=1
            string = string + new_string

        elif i < len(lista)-1 :
            if len(lista[0:npop-pop]) >=2:
                new_string = ','+str(pop)+':'+str(sum(lista[thesi:thesi+len(lista[0:npop-pop])])/len(lista[0:npop-pop]))
                thesi= thesi +len(lista[0:npop-pop])
                pop+=1
                string = string + new_string
            
            
        elif i == len(lista)-1 or i == len(lista):
            new_string = ','+str(pop)+':'+str(lista[-1])
            string = string + new_string
            i+=1
            pop+=1
            

    string = '('+string+');'
    return string      


def RF_tree_distance(old_time_string,after_pca_time_string):
    ## Finds RF distance between two trees

    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(
            data=old_time_string,
            schema='newick',
            taxon_namespace=tns)
    tree2 = dendropy.Tree.get(
            data=after_pca_time_string,
            schema='newick',
            taxon_namespace=tns)
    RF = (treecompare.weighted_robinson_foulds_distance(tree1, tree2),'rf') #check symmetry of trees and based on lentghs
    
    return RF


def fst_similarity(fst_dist,each_old_distance,number_of_pop):
    ## Compares fst index with scenario population distance
    
    sumofel_before = 0
    sumofel_after = 0
    elements1 = 0
    elements2 = 0
    anw_before = each_old_distance
    anw_before = each_old_distance[np.triu_indices(3, k = 1)]
    anw_after = np.asarray(fst_dist)
    for i in anw_before :
        if i!=min(anw_before):
            sumofel_before = sumofel_before +i
        if i==min(anw_before):
            elements1+=1

    if elements1 == 2:
        print('there are elements with same value')
        return 0 
    else:
        mean_of_anw_before = sumofel_before/2 #find mean of 2 biggest elements

    for i in anw_after :
        if i!=min(anw_after):
            sumofel_after = sumofel_after +i
        if i==min(anw_after):
            elements2+=1

    if elements2 == 2 or elements2==3: ##CHANGED
        print('there are elements with same value')
        return 0 
    else:
        mean_of_anw_after = sumofel_after/2 #find mean of 2 biggest elements

    minof_before = min(anw_before)
    trimmed_list_before = [x for x in anw_before if x!=minof_before]  
    diafora_trimmed_before = abs(trimmed_list_before[0]-trimmed_list_before[1])
    diafora_small_numbers_before = abs(min(trimmed_list_before)-minof_before)

    if diafora_small_numbers_before < diafora_trimmed_before:
        return 0 

    minof_after = min(anw_after)
    trimmed_list_after = [x for x in anw_after if x!=minof_after]  
    diafora_trimmed_after = abs(trimmed_list_after[0]-trimmed_list_after[1])
    diafora_small_numbers_after = abs(min(trimmed_list_after)-minof_after)
    
    if diafora_small_numbers_after < diafora_trimmed_after:
        return 0
        
    if (mean_of_anw_before > min(anw_before)) and  (mean_of_anw_after > min(anw_after)) :
        if np.argmin(anw_before) == np.argmin(anw_after):
            return 1
        else:
            return 0 
    else:
        return 0 


def new_fst(genot_matrix,fst_after_reps,number_of_pop):
    """ Calculates nucleotide divergence
    
    Parameters
    ----------
    genot_matrix: arraylike, matrix with genotypes rows:snps ,columns: number of samples
    fst_after_reps: list, calculated fst index
    number_of_pop: int, number of populations

    Returns
    ------
    fst_after_reps: list, contains all fst indexes after all reps
    fst_between_pops: list, contains fst indexes in each rep
    """
    merge_genome = []
    fst_between_pops = []
    
    print(genot_matrix.shape,' SHAPE_DATA_AFTER_FST')
    genot_matrix = genot_matrix.T #Samples x Snps
    GENOME_MATRIX_ALL_POPS = np.split(genot_matrix,number_of_pop)
    
    for a, b in itertools.combinations(GENOME_MATRIX_ALL_POPS, 2): #
        a = a.T 
        b = b.T
        
        merge_pop_genome = np.concatenate((a, b), axis=1) #snps x samples
        merge_genome.append(merge_pop_genome)
        
    for genome in merge_genome:
        #print(type(genome),genome,genome.shape,'genome\n')
        _,columns = genome.shape
        haplotypes = columns/2
        genome = allel.HaplotypeArray(genome)
        #print('genome',genome,'\n')
        ac1 = genome.count_alleles(subpop=[0,haplotypes-1])
        ac2 = genome.count_alleles(subpop=[haplotypes,columns-1])
        num, den = allel.hudson_fst(ac1, ac2)
        fst = np.sum(num) / np.sum(den)
        fst_between_pops.append(fst)
    fst_after_reps.append(fst_between_pops)

    return fst_after_reps,fst_between_pops

def plot_parameters_scores(param_list,samples_values,supplement_param,config):
    #print(param_list,'MPAMP')
    
    """ Finds which parameters are about to be plotted
    
    Parameters
    ----------
    param_list: list, names of variable parameters
    samples_values: list, values if samples of each population
    supplement_param: SupplementaryParameter, contains supplementary parameters of simulation
    config: dict, contains all parameters of simulation

    Returns
    ------
    MATRIX: arralike, ld_clumped mat with rows:snps, columns:number_of_samples
    """
    param_combination_list = []
    if len(param_list) == 1:
        if param_list[0] == 'sample':
            return samples_values
        elif param_list[0] == 'pop_config':
            return supplement_param.population_configurations
        elif param_list[0] == 'migration':
            return config['mig']
        elif param_list[0] == 'mut_rate':
            return config['mutation_rate']
        elif param_list[0] == 'length':
            return config['length']
        elif param_list[0] == 'recomb_rate':
            return config['recombination_rate']
        elif param_list[0] == 'demo_event':
            return supplement_param.demographic_events
        elif param_list[0] == 'distance':
            return config['distance']
        elif param_list[0] == 'bottleneck':
            return supplement_param.bottle_end_list[0]
        elif param_list[0] == 'recent_merge':
            return supplement_param.recent_merge_list[0]
        elif param_list[0] == 'oldest_merge':
            return supplement_param.oldest_merge_list[0]

    elif len(param_list)==2:
        if  'sample' in param_list:
            param_combination_list.append( samples_values)
        if 'pop_config' in param_list:
            param_combination_list.append( supplement_param.population_configurations)
        if 'migration' in param_list:
            param_combination_list.append(config['mig'])
        if 'mut_rate' in param_list:
            param_combination_list.append(config['mutation_rate']) 
        if 'length' in param_list:
            param_combination_list.append(config['length']) 
        if 'recomb_rate' in param_list:
            param_combination_list.append(config['recombination_rate'])  
        if 'demo_event' in param_list:
            param_combination_list.append(supplement_param.demographic_events)   
        if 'distance' in param_list:
            param_combination_list.append(config['distance'])
        if 'bottleneck' in param_list:
            param_combination_list.append([supplement_param.bottle_end_list[0]])
        if 'recent_merge' in param_list:
            
            
            param_combination_list.append([supplement_param.recent_merge_list[0]])
        if 'oldest_merge' in param_list:
            param_combination_list.append([supplement_param.oldest_merge_list[0]])
        
        return param_combination_list
    else:
        print('no valid parameters given')
        print(param_list,len(param_list))
        exit()

def ld_clumping(mat):
    """ Reduces number of SNPs by conducting LD clumping 
    
    Parameters
    ----------
    mat: arraylike, rows:snps, columns:number_of_samples

    Returns
    ------
    MATRIX: arralike, ld_clumped mat with rows:snps, columns:number_of_samples
    """
    MAF_LIST = []
    mat_transpose = mat
    
    for row,i in enumerate(mat_transpose):
        if 2 not in mat:
            count_of_1 = np.count_nonzero(i == 1)/len(i)
            count_of_0 = np.count_nonzero(i == 0)/len(i)
            
            if count_of_0 < count_of_1:
                MAF = count_of_0 #minimum allel frequency
                MAF_LIST.append((MAF,row))
            else:
                MAF = count_of_1
                MAF_LIST.append((MAF,row))
        elif 2 in mat:
            count_of_2 = np.count_nonzero(i == 2)/len(i)
            count_of_1 = np.count_nonzero(i == 1)/len(i)
            count_of_0 = np.count_nonzero(i == 0)/len(i)
            
            if count_of_0 < count_of_1 or count_of_0 < count_of_2:
                MAF = count_of_0 #minimum allel frequency
                MAF_LIST.append((MAF,row))
            elif count_of_1 < count_of_0 or count_of_1 < count_of_2:
                MAF = count_of_1
                MAF_LIST.append((MAF,row))
            else:
                MAF = count_of_2
                MAF_LIST.append((MAF,row))

    diafora = 1
    if diafora ==1:
        indexes_of_sorted_MAF = [j[1] for j in sorted(MAF_LIST ,key=itemgetter(0), reverse=True)] #MAF_LIST
        diafora=2
   
    cur_index = 0
    while cur_index < len(indexes_of_sorted_MAF)-1:
        first_snp = indexes_of_sorted_MAF[cur_index] #index of first snp in index sorted maf list,having the biggest MAF
        second_snp = indexes_of_sorted_MAF[cur_index+1] #index of second snp in index sorted maf list,having the biggest MAF
        a = list(mat_transpose[first_snp]) #snp1
        b = list(mat_transpose[second_snp]) #snp2
        
        cor = np.corrcoef(a,b)
        r = cor[0,1]
        
        if r >= 0.2 : #0.5
            #diafora+=1
            del indexes_of_sorted_MAF[cur_index+1]
        else:
            cur_index = cur_index+1

    sorted_indexes = sorted(indexes_of_sorted_MAF)
    MATRIX =  mat_transpose[sorted_indexes, :]
    
    print('shape of final after ld_clumping',MATRIX.shape)
    return MATRIX


def writing_scores(supplement_param,scor_di,parameter,parameter_values,folder,samples_values): #writes fst,euclideian and nucl_diff scores in output files
    """ Preprocessing of data that are about to be written in ouptut files
    
    Parameters
    ----------
    supplement_param: SupplementaryParameter, contains supplementary parameters of simulation
    scor_di: dict, contains lists of different type of scores from recognition process of populations
    parameter: list, contains names of variables of simulation
    parameter_values: list, contains values of variables of simulation
    folder: str, the directory where all files are kept
    samples_values: list, values if samples of each population
    """
    #print(parameter_values,'breeeeeegiatidenexei paarmeters')
    migr_pops_matrix = []
    migr_pops_list = []
    migr_ratio_list = list()
    combine_values = []
    mig_li = []
    
    s = [str(i+1) for i in supplement_param.pop_bottled_list]
    
    
    supplement_param.pop_bottled_list = ",".join(s)
    
    
    path_of_scores_to_plot = folder+'scores_to_plot'
    createFolder(path_of_scores_to_plot)
    
    if type(parameter) == str:
        parameter = [parameter]
    
    if len(parameter) == 1:
        migration_population_list = []
        migration_pop = ''
        
        if parameter[0] == 'migration':
            for mig_matrix in  parameter_values:
                mig_matrix = np.asarray(mig_matrix)
                migr_pops,migr_ratio = find_migration_ratios(mig_matrix)
                migr_pops_list.append(migr_pops)
                migr_ratio_list.append(migr_ratio)
    
            parameter_values = migr_ratio_list
            for i in migr_pops_list:
                if i == [0]:
                    migration_population_list.append('0')
                else:
                    for every_migr_pop in i :
                        if  migration_pop == '':
                            migration_pop = (str(every_migr_pop[0])+','+str(every_migr_pop[1]))
                            mig_li.append(str(every_migr_pop[0])+','+str(every_migr_pop[1]))
                        else:
                            if (str(every_migr_pop[0])+','+str(every_migr_pop[1])) not in mig_li:
                                migration_pop = migration_pop +'|'+ (str(every_migr_pop[0])+','+str(every_migr_pop[1]))
                                mig_li.append((str(every_migr_pop[0])+','+str(every_migr_pop[1])))
                    migration_population_list.append(migration_pop)
            migr_pops_matrix = migration_population_list
            
    elif len(parameter) == 2 :
        if 'migration' in parameter:
            #print('MI1')
            for i,param in enumerate(parameter):
               
                if param == 'migration':
                    
                    for mig_matrix in  parameter_values[i]:
                        mig_matrix = np.asarray(mig_matrix)
                        migr_pops,migr_ratio = find_migration_ratios(mig_matrix)
                        migr_pops_list.append((migr_pops))
                        migr_ratio_list.append(migr_ratio)
                        
                    parameter_values[i] = migr_ratio_list
                # else: 
                    repeat = len(parameter_values[i])
                    migration_population_list = []
                    migration_pop = ''
                    migration_pop_set = set()
                    for i in migr_pops_list:
                        if i == [0]:
                            migration_population_list.append('0')
                        else:
                            for every_migr_pop in i :
                                if  migration_pop == '':
                                    mig_li.append(str(every_migr_pop[0])+','+str(every_migr_pop[1]))
                                    migration_pop = (str(every_migr_pop[0])+','+str(every_migr_pop[1]))
                                else:
                                    if (str(every_migr_pop[0])+','+str(every_migr_pop[1])) not in mig_li:
                                        #if migration_pop not in migration_pop_set:
                                        migration_pop_set.add(migration_pop)
                                        migration_pop = migration_pop +'|'+ (str(every_migr_pop[0])+','+str(every_migr_pop[1]))
                                        mig_li.append((str(every_migr_pop[0])+','+str(every_migr_pop[1])))
                            migration_population_list.append(migration_pop)    
                    
                    if  ('mut_rate' in parameter) or ('length' in parameter) or ('recomb_rate' in parameter):
                        migr_pops_matrix = np.repeat(migration_population_list,repeat)
                        
                    else:
                        migr_pops_matrix = migration_population_list*repeat
        
        
        
        if event == 'bottle_recent_merge' or event == 'bottle_oldest_merge':
            print('mpika logw tou bottlerecentmerge\n')
            combine_values = [(i,j) for i in parameter_values[0] for j in parameter_values[1]]
        else:
            combine_values = [(i,j) for i in parameter_values[0] for j in parameter_values[1] if i != j]
       
        
    for number,param in enumerate(parameter):
        if len(parameter) == 2 and param == 'sample':
            parameter_values[number] = samples_values
        else:
            continue
    
        combine_values = [(i,j) for i in parameter_values[0] for j in parameter_values[1] if i != j]
        
    
    path = folder+'scores_to_plot/'
    scor_par = ScoresParam()
    scor_par.parameter = parameter
    scor_par.parameter_values = parameter_values
    scor_par.path = path
    scor_par.migr_pops_matrix = migr_pops_matrix
    scor_par.combine_values = combine_values
    
    for key,_ in scor_di.items():
        name_of_out_file = str(key)+'diff_success.txt'
        
        write_scores_to_output_file(name_of_out_file,scor_par,supplement_param,scor_di[key])
    

    
def computation_of_scores(scores_lists,scor_di):
    """ Modifies scores in *100 unit
    
    Parameters
    ----------
    scores_lists: dict, contains lists of number of success (0,1) from recognition process of populations
    scor_di: dict, contains lists of different type of scores from recognition process of populations

    Returns
    ------
    scor_di
    """
    
    for scr_li in scores_lists.keys():
        suc = scores_lists[scr_li].count(1)/len(scores_lists[scr_li])
        suc_portion = int(suc *100)
        scor_di[scr_li].append(suc_portion)

    return scor_di

def find_migration_ratios(migr_matrix):
    """ Seeks for migration ratios from input data
    
    Parameters
    ----------
    migr_matrix: arraylike, where all migration ratios are

    Returns
    ------
    migr_pops: tuple, numbers of populations
    migr_ratio[0]: float, migration ratio of migration between populations in simulation
    """

    all_zeros = np.nonzero(np.any(migr_matrix != 0, axis=0))[0]
    
    migr_pops = []
    if len(all_zeros) == 0:
        
        migr_ratio = 0
        migr_pops.append((0))
        return migr_pops,migr_ratio
    if len(all_zeros) != 0:
        zero = np.nonzero(migr_matrix)
        migr_ratio = migr_matrix[zero]
        for i in zip(zero[0], zero[1]):
            pop_a = int(i[0])+1
            pop_b = int(i[1])+1
            if ((pop_b,pop_a)) not in migr_pops:
                migr_pops.append((pop_a,pop_b))

        return migr_pops,migr_ratio[0]
    

def print_header(f,parameter,parameter_values):
    """ Prints header in output files
    
    Parameters
    ----------
    f: file, where output is about to be written
    parameter: list, contains names of variables of simulation
    parameter_values: list, contains values of variables of simulation
    """

    if len(parameter)==1:
        if parameter[0] != 'bottleneck':
            if parameter[0]== 'recent_merge' or parameter[0] == 'oldest_merge':
                f.write('{} merged_pops succ_score\n'.format(parameter[0]))
            elif parameter[0]=='migration':
                f.write('{} migr_pops succ_score\n'.format(parameter[0]))
            else:
                f.write('{} succ_score\n'.format(parameter[0]))
        else:
            f.write('t_start t_end_of bottleneck_pops succ_score\n')

    if len(parameter)==2:
        if 'bottleneck' in parameter and 'migration'in parameter:
            if parameter[0] == 'bottleneck':
                f.write('{} {} {} {} {} {} \n'.format('bottle_end',parameter[1],'bottle_st','pop_bottled','migr_pops','succ_score'))
            else:   
                f.write('{} {} {} {} {} {} \n'.format(parameter[0],'bottle_end','bottle_st','pop_bottled','migr_pops','succ_score'))
            
        elif ('bottleneck' in parameter and 'recent_merge'in parameter) or ('bottleneck' in parameter and 'oldest_merge'in parameter):
            if parameter[0] == 'bottleneck':
                f.write('{} {} {} {} {} {} \n'.format('bottle_end',parameter[1],'bottle_st','pop_bottled','merged_pops','succ_score'))
            else:
                f.write('{} {} {} {} {} {} \n'.format(parameter[0],'bottle_end','bottle_st','pop_bottled','merged_pops','succ_score'))
                
        elif ('migration' in parameter) and ('recent_merge'in parameter) or ('migration' in parameter) and ('oldest_merge'in parameter):
            f.write('{} {} {} {} {} \n'.format(parameter[0],parameter[1],'migr_pops','merged_pops','succ_score'))
        else: 
            if ('bottleneck' not in parameter) and ('recent_merge' not in parameter) and ('oldest_merge' not in parameter) and ('migration'not in parameter):
                f.write('{} {} {} \n'.format(parameter[0],parameter[1],'succ_score'))
            elif ('bottleneck' in parameter) and ('recent_merge' not in parameter) or ('oldest_merge' not in parameter) and ('migration' not in parameter):
                if parameter[0] == 'bottleneck':
                    f.write('{} {} {} {} {} \n'.format('bottle_end',parameter[1],'bottle_st','pop_bottled','succ_score'))
                else:
                    f.write('{} {} {} {} {} \n'.format(parameter[0],'bottle_end','bottle_st','pop_bottled','succ_score'))

            elif ('recent_merge' in parameter) or ('oldest_merge' in parameter) and ('bottleneck' not in parameter) and ('migration'not in parameter):
                f.write('{} {} {} {} \n'.format(parameter[0],parameter[1],'merged_pops','succ_score'))
                
            elif ('migration'in parameter) and ('bottleneck' not in parameter) and ('recent_merge' not in parameter) and ( 'oldest_merge' not in parameter):
                f.write('{} {} {} {} \n'.format(parameter[0],parameter[1],'migr_pops','succ_score'))


            

            
def hudson_fst(ac1, ac2, fill=np.nan):
    """Calculate the numerator and denominator for Fst estimation using the
    method of Hudson (1992) elaborated by Bhatia et al. (2013).

    Parameters
    ----------
    ac1 : array_like, int, shape (n_variants, n_alleles)
        Allele counts array from the first population.
    ac2 : array_like, int, shape (n_variants, n_alleles)
        Allele counts array from the second population.
    fill : float
        Use this value where there are no pairs to compare (e.g.,
        all allele calls are missing).

    Returns
    -------
    num : ndarray, float, shape (n_variants,)
        Divergence between the two populations minus average
        of diversity within each population.
    den : ndarray, float, shape (n_variants,)
        Divergence between the two populations.

    """

    # # check inputs
    ac1 = asarray_ndim(ac1, 2)
    ac2 = asarray_ndim(ac2, 2)
    check_dim0_aligned(ac1, ac2)
    ac1, ac2 = ensure_dim1_aligned(ac1, ac2)

    # calculate these once only
    an1 = np.sum(ac1, axis=1)
    an2 = np.sum(ac2, axis=1)

    # calculate average diversity (a.k.a. heterozygosity) within each
    # population
    within = (mean_pairwise_difference(ac1, an1, fill=fill) +
              mean_pairwise_difference(ac2, an2, fill=fill)) / 2

    # calculate divergence (a.k.a. heterozygosity) between each population
    between = mean_pairwise_difference_between(ac1, ac2, an1, an2, fill=fill)

    # define numerator and denominator for Fst calculations
    num = between - within
    den = between

    return num, den

def subtract_float(a,b):
    """ Substracts numbers 
    
    Parameters
    ----------
    a: int/float number
    b: int/float number
    
    Returns
    -------
    difference: Difference between numbers
    """
    
    x = Decimal(str(a))
    y = Decimal(str(b))
    difference = float(x-y)
    return difference

def multiple_pops_similarity(each_distance_before,distances_of_pops,number_of_pop):
    """ Compares old distance matrix to the new one(after pca,fst etc) between multiple populations
    
    Parameters
    ----------
    each_distance_before: arraylike, old distance matrix of populations (scenario distance)
    distances_of_pops: arraylike, new distance matrix (after pca,fst,tsne process)
    number_of_pop: int ,the number of compared populations 
    
    Returns
    -------
    0/1: 0 if old distance matrix does not have the same structure with new distance matrix. 1 if there is same structure, 
    so populations can be recognised 
    """
    each_distance_before = each_distance_before[np.triu_indices(number_of_pop, k = 1)]
    distances_of_pops = np.asarray(distances_of_pops)
    tri = np.zeros((number_of_pop,number_of_pop))
    tri[np.triu_indices(number_of_pop,1)] = distances_of_pops
    

    before_result_max = np.where(each_distance_before == np.amax(each_distance_before))
    before_result_min = np.where(each_distance_before == np.amin(each_distance_before))
    after_result_min = np.where(distances_of_pops == np.amin(distances_of_pops))
    
    
    max_elements = []
    
    if len(before_result_min[0])==1 and len(after_result_min[0])==1 and before_result_min == after_result_min:
        for i in before_result_max[0]:
            max_elements.append(distances_of_pops[i])
        
        for pos,row in enumerate(tri[:-1]):
            max_difference = -1
            m = np.min(row[np.nonzero(row)])
        
            for i in row:
                diff=subtract_float(i,m)
                if diff > max_difference:
                    max_difference = diff
    
            for next_row in tri[pos+1:]:
                max_next = max(next_row)
                if max_next > m: #if max(next_row) > min(previous_row)
                    return 0 
                diff=abs(subtract_float(max_next,m))
                
                
                if diff <= max_difference:
                    return 0
        return 1
    else:
        print('eixe polles mikres times')
        return 0






def print_new_header(file,parameter):
    """ Prints header in output files
    
    Parameters
    ----------
    f: file, where output is about to be written
    parameter: list, contains names of variables of simulation
    """
    
    with open(file, 'a') as f:
        if len(parameter)==1:
            p = []
            p = [parameter[0],'merged_pops','succ_score']
            if parameter[0] != 'bottleneck':
                if parameter[0]== 'recent_merge' or parameter[0] == 'oldest_merge':
                    pass
                    
                elif parameter[0]=='migration':
                    
                    p[1] = 'migr_pops'
                else:
                    
                    p[1] = 'succ_score'
                    p[2] = ''
                    
            else:
                p = ['bottle_st','bottle_end','pop_bottled','succ_score']
                

        if len(parameter)==2:
            p=[]
            if 'bottleneck' in parameter:
                p = [parameter[0],'bottle_end','bottle_st','pop_bottled','migr_pops','succ_score']
                if parameter[0] == 'bottleneck':
                    p[0] = "bottle_end"
                    p[1] = parameter[1]
                if 'recent_merge'in parameter or 'oldest_merge'in parameter:
                    p[4] = 'merged_pops'
                elif ('recent_merge' not in parameter) and ('oldest_merge' not in parameter) and ('migration' not in parameter):
                        if parameter[0] == 'bottleneck':
                            p[0] = 'bottle_end'
                            p[1] = parameter[1]
                        else:
                            p[1] = 'bottle_end'
                        p[4] = 'succ_score'
                        p[5] = ''
            else:
                p = [parameter[0],parameter[1]]
                if ('migration' in parameter):
                    if ('recent_merge'in parameter) or ('oldest_merge'in parameter):
                        p.extend(['migr_pops','merged_pops','succ_score'])
                    elif ('recent_merge' not in parameter) and ( 'oldest_merge' not in parameter):
                        p.extend(['migr_pops','succ_score'])
                else:
                    if ('recent_merge' not in parameter) and ('oldest_merge' not in parameter): 
                        p.append('succ_score')

                    elif (('recent_merge' in parameter) or ('oldest_merge' in parameter)): 
                        p.extend(['merged_pops','succ_score'])


        for i in range(len(p)):
            if i=='':
                pass
            else:
                f.write("{} ".format(p[i]))
        f.write("\n")
    f.close()
    
    

    return p


def write_scores_to_output_file(output_name,scor_par,supplement_param,SUCC_SCORES):
    path = scor_par.path
    filename = path+output_name

    already_exists = os.path.exists(filename)

    
    data = {}
    
    with open(filename,'a+') as myfile:
        if not already_exists:
            
            header = print_new_header(filename,scor_par.parameter)
        
        if already_exists:
            with open(filename,'r') as f:
                firstline = f.readline().rstrip('\n').split(' ')[:-1]
                header = firstline
            f.close
            
        if type(scor_par.parameter_values) != list:
            scor_par.parameter_values = [scor_par.parameter_values]

        if type(scor_par.parameter_values) == list:

            if type(supplement_param.merged_pops) == tuple:
                supplement_param.merged_pops = str(supplement_param.merged_pops[0])+','+str(supplement_param.merged_pops[1])
            data["bottle_st"] = supplement_param.bottle_start
            data["pop_bottled"] = supplement_param.pop_bottled_list
            data["merged_pops"] = supplement_param.merged_pops
            
            loop_var = scor_par.parameter_values
            
            if len(scor_par.parameter)==2:
                loop_var = scor_par.combine_values
            
            for num in range(0,len(loop_var)): #parameters
                data["succ_score"]=SUCC_SCORES[num]
                values_to_print = loop_var
                
                if len(scor_par.parameter)==2:
                    values_to_print = loop_var[num]
                    
                    
                if len(scor_par.parameter)==1:
                    if "bottle_st" in header or scor_par.parameter[0] =='bottleneck':
                        val_key = "bottle_end"
                    elif "migration" in header:
                        
                        data["migr_pops"] = scor_par.migr_pops_matrix[num]
                        
                        val_key = "migration"
                    elif "recent_merge" in header:
                        val_key = "recent_merge"
                    elif "oldest_merge" in header:
                        val_key = "oldest_merge"
                    else:
                        val_key = header[0]
                    data[val_key] = values_to_print[num]
                if len(scor_par.parameter)==2:
                    if 'bottleneck' in scor_par.parameter:
                        if 'migration'in scor_par.parameter:
                            val_key2 = "bottle_end"
                            val_key1 = "migration"
                            
                            data["migr_pops"] = scor_par.migr_pops_matrix[num]
                        elif ('recent_merge'in scor_par.parameter)  or ('oldest_merge'in scor_par.parameter):
                            val_key1 = "bottle_end"
                            val_key2 = "oldest_merge"
                            if "recent_merge" in header:
                                val_key2 = "recent_merge"
                            
                        elif ('recent_merge' not in scor_par.parameter) and ('oldest_merge' not in scor_par.parameter) and ('migration' not in scor_par.parameter):
                            val_key1 = header[0]
                            val_key2 = header[1]
                    else:
                        if 'migration' in scor_par.parameter:
                            if ('recent_merge'in scor_par.parameter) or ('oldest_merge'in scor_par.parameter):
                                data["migr_pops"] = scor_par.migr_pops_matrix[num]
                                val_key1 = "migration"
                                val_key2 = "oldest_merge"
                                if "recent_merge" in header:
                                    val_key2 = "recent_merge"
                                
                            elif ('recent_merge' not in scor_par.parameter) and ( 'oldest_merge' not in scor_par.parameter):
                                data["migr_pops"] = scor_par.migr_pops_matrix[num]
                                val_key1 = header[0]
                                val_key2 = header[1]
                        else:
                            if ('recent_merge' not in scor_par.parameter) and ( 'oldest_merge' not in scor_par.parameter):
                                
                                val_key1 = header[0]
                                val_key2 = header[1]
                
        
                            elif ('recent_merge' in scor_par.parameter) or ('oldest_merge' in scor_par.parameter):
                                val_key1 = header[0]
                                val_key2 = header[1]
                                
                    data[val_key1] = values_to_print[0]
                    data[val_key2] = values_to_print[1]

                for column in header:
                    if column == '':
                        pass
                    else:
                        myfile.write('{} '.format(data[column]))
                myfile.write('\n')

    myfile.close()

def make_distance(dist_el):
    maximum = 0
    for i in dist_el:
        numbers = i.split('=')[0].split(',')
        numbers = [int(x) for x in numbers]
        max_num = max(numbers)
        if max_num > maximum:
            maximum = max_num
        
    number_of_pops = maximum
    Distance = np.zeros((number_of_pops,number_of_pops))



    for i in dist_el:
        pops = i.split('=')[0].split(',')
        dist = i.split('=')[1]
        
        Distance[int(pops[0])-1][int(pops[1])-1] = int(dist)
        
    
    return Distance,number_of_pops




def make_migration_matrix(mig,number_of_pops):
    mig_mat_list = []
    if len(mig) ==1 and mig[0]=='':
        migration_matrix = np.zeros((number_of_pops,number_of_pops))
        migr_mat_list = [[ matrix.tolist() for matrix in migration_matrix]]
        
        return migr_mat_list
    else:
        for i in mig:
            if i =='':
                migration_matrix = np.zeros((number_of_pops,number_of_pops))
                mig_mat_list.append(migration_matrix)
            else:
                migpops = i.split('=')[0].split('/')
                migration_matrix = np.zeros((number_of_pops,number_of_pops))
                for mi in migpops:
                    mig_rate = i.split('=')[1]
                    migration_matrix[int(mi[0])][int(mi[2])] = float(mig_rate)
                mig_mat_list.append(migration_matrix)
                
                

        return mig_mat_list










if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script uses a txt file with parameters as input to simulate via msprime a scenario about a number of populations. It uses PCA in order to reduce dimensions of the msprime produced data and calculates distances between populations per simulation. Afterwards there is a check of similarity between populations distance before and after PCA')
    parser.add_argument('prefix', help='prefix for output folder')
    
    parser.add_argument('-event', help='bottleneck effect,merge of pops') 
    parser.add_argument('-DE', help='demographic events',action='append') 
    parser.add_argument('-samples', help='samples of pops',action='append') 
    parser.add_argument('-N', help='effective population size') 
    parser.add_argument('-R', help='population growth rate') 
    parser.add_argument('-length', help='length of sequence',action='append') 
    parser.add_argument('-mut', help='mutation rate',action='append') 
    parser.add_argument('-recomb', help='recombination rate',action='append') 
    parser.add_argument('-mig', help='migration rate',action='append') 
    parser.add_argument('-dist', help='distances of populations in years',action='append') 


    


    
    
    args = parser.parse_args()
    prefix = args.prefix
    
    event = args.event
    DE = args.DE
    simvar.samples = args.samples
    simvar.N = [int(n) for x in [args.N]  for n in x.split(' ')]
    simvar.R = [float(r) for x in [args.R]  for r in x.split(' ')]
    simvar.length = args.length
    simvar.mutation_rate = args.mut
    simvar.recombination_rate = args.recomb

    simvar.mig = args.mig
    simvar.dist = args.dist
    

    
    main()
    print('Total time: ',datetime.now() - startTime)
    
