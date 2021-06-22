import os
import argparse
import glob
import sys
import itertools
from datetime import datetime



def multiply_param(times, value):
    return (times * value)[:-1]

# From command line
classic_merge = 'empty'
param1 =  sys.argv[1]#"recent"
param2 = sys.argv[2]#"samples"
number_of_pops = int(sys.argv[3])
timestring = sys.argv[4].split(',') #string with time of merges
event_factor = sys.argv[5]
migration_factor = sys.argv[6]
bottle_factor = sys.argv[7]     


variables = {}
variables["length"] = [2e5,4e5,6e5,8e5  ,10e5]
variables["samples"] = [20,40,100,160,200]
variables["recomb"] = [2e-9,2e-8,2e-7,2e-6]
variables["mut"] = [1e-9,1e-8,1e-7,1e-6]#1e-8,1e-7,1e-6,
variables["N"] = multiply_param(number_of_pops, '10000 ')
variables["R"] = multiply_param(number_of_pops, '0 ')
variables["recent_merge"] = ['']
variables["oldest_merge"] = ['']
variables["bottleneck"] = ['']
variables["migration"] = ['']

default = {}
default["length"] = 1e5
default["samples"] = 100
default["recomb"] = 2e-8
default["mut"] = 1e-8
default["N"] = multiply_param(number_of_pops, '10000 ')
default["R"] = multiply_param(number_of_pops, '0 ')
default["recent_merge"] = ''
default["oldest_merge"] = ''
default["bottleneck"] = ''
default["migration"] = ''


def arg_string(key, value):
    return '-{} "{}" '.format(key, value)


def get_pairs(number_of_pops):
    return list(itertools.combinations(range(1,number_of_pops+1),2))

def is_merge(i, j):
    return (i=="recent") or (i=="oldest") or (j=="recent") or (j=="oldest")


def new_arg_merge_string(number_of_pops,distances):
    
    num=0
    dists = ""
    for pair in get_pairs(number_of_pops):
        
        dists = dists + ' -dist "{},{}={}"'.format(pair[0],pair[1],distances[num])
        num+=1
    return dists


def check_params(param1,param2):
    if (param1 =='recent_merge' or param1 =='oldest_merge' or param1 =='bottleneck' or param1 =='migration') and (param2 =='recent_merge' or param2 =='oldest_merge' or param2 =='bottleneck' or param2 =='migration'):
        return True

def complicate_not_exist(var):
    if var =='recent_merge' or var =='oldest_merge' or var =='bottleneck' or var =='migration':
        return 'complicated'
def complicate_par(bottlenecks,migration,merge_t):   
    for i in migration:
        migr= i.split(' ')
        migr_count = migr.count('-mig')
        

    if bottlenecks!=[''] or migr_count>1 or merge_t!=['']:
        
        return True
    else:
        return False





def bottleneck(bottle_factor):
    
    if bottle_factor !='':
        not_range_bottlenecks = []
        bottlenecks = []
        fact_string = ""
        time_change = 0
        size_change = 0
        size_change_facts = bottle_factor.split('/')
        for s in size_change_facts:
            fact = s.split(":")
            if '-' in fact[1]: #range in time
                
                time_range = []
                start = int(fact[1].split('-')[0])
                end = int(fact[1].split('-')[1])
                step = int(fact[1].split('-')[2])
                time_change = 1
                if ',' in fact[3]:
                    populations = fact[3].split(',')
                    
                    for time in range(start,end+step,step):
                        middle_str = ''
                        for p in populations:
                            middle_str = middle_str+ '-DE "de {} {} {} {}" '.format(fact[0],time,fact[2],int(p))
                        time_range.append(middle_str)
                    
                else:
                    for time in range(start,end+step,step):
                        time_range.append(' -DE "de {} {} {} {}" '.format(fact[0],time,fact[2],fact[3]))
                    
            if '-' in fact[2]: #range in pop size
                size_range = []
                st_size =  int(fact[2].split('-')[0])
                end_size =  int(fact[2].split('-')[1])
                step_size = int(fact[2].split('-')[2])
    
                size_change = 1
                if ',' in fact[3]:
                    populations = fact[3].split(',')
                    for size in range(st_size,end_size+step_size,step_size):
                        middle_size_str = ''
                        for p in populations:
                            middle_size_str = middle_size_str+ ' -DE "{} {} {} {}" '.format(fact[0],fact[1],size,int(p))
                        size_range.append(middle_size_str)
                else:
                    for size in range(st_size,end_size+step_size,step_size):
                        size_range.append(' -DE "de {} {} {} {}" '.format(fact[0],fact[1],size,fact[3]))
            
            if '-' not in fact[1] and '-' not in fact[2]:
                if ',' in fact[3]:
                    populations = fact[3].split(',')
                    for po in populations:
                        fact_string = fact_string+ ' -DE "de {} {} {} {}" '.format(fact[0],fact[1],fact[2],int(po))
                    not_range_bottlenecks.append(fact_string)
                else:
                    fact_string = fact_string+ ' -DE "de {} {} {} {}" '.format(fact[0],fact[1],fact[2],fact[3])
                    not_range_bottlenecks.append(fact_string)

        if time_change != 0:
            for each_time in time_range:
                for each_bottle_event in not_range_bottlenecks:
                    bottlenecks.append(each_bottle_event+each_time)
        elif size_change != 0:
            for each_size in size_range:
                for each_bottle_event in not_range_bottlenecks:
                    bottlenecks.append(each_bottle_event+each_size)
        else:
            for each_bottle_event in not_range_bottlenecks:
                    bottlenecks.append(each_bottle_event)


    else:
        bottlenecks = [''] 
    return bottlenecks

def mass_migration(timestring,pops,param1,param2): 
    
    num = 0
    time_merge_list = []
    mass_pops = []
    final_mass = []

    
    if (param1 != "oldest_merge") or (param1 == "recent_merge") and (param2 == "oldest_merge") or (param2 != "recent_merge"):
        for k in range(pops-1,0,-1):
            l = k*(str(timestring[num])+' ')
            time = l.split(' ') #
            if len(time_merge_list) >1:
                
                time_merge_list = time_merge_list+[t for t in time if t !='']
            else:
                time_merge_list = [t for t in time if t !='']
            num+=1
            time_merge_list = [int(y) for y in time_merge_list]
        classic_merge = new_arg_merge_string(pops,time_merge_list)
        #print(classic_merge,'classic_merge')
    else:
        classic_merge = 'None'
        


    timestring = [int(x) for x in set(timestring)]
    sort_timestring = sorted(timestring)
    
    for i in range(1,pops+1):
        if i !=1:
            mass_pops.append((i,i-1))
    mass_pops = sorted(mass_pops,reverse=True)
    #print(mass_pops,'MASS')
    
    mass_mig_string = ""
    for n,_ in enumerate(sort_timestring):
        mass_mig_string = mass_mig_string + ' -DE "de mm {} {} {} 1.0"'.format(sort_timestring[n],int(mass_pops[n][0])-1,int(mass_pops[n][1])-1)
    final_mass.append(mass_mig_string)
    return final_mass,classic_merge
    


def process_merge(timestring,pops,merge_range,ti,time_merge_list,num,type_merge,mass_events_list):
    for i in range(pops-1,0,-1):
        
        if type_merge == "oldest":
            timestring[0] = str(ti)
        elif type_merge =='recent':
            timestring[-1] = str(ti)

        l = i*(timestring[num]+' ')
        time = l.split(' ') 
        if len(time_merge_list) >1:
            
            time_merge_list = time_merge_list+[t for t in time if t !='']
        else:
            time_merge_list = [t for t in time if t !='']
        num+=1
    time_merge_list = [int(y) for y in time_merge_list]

    mass_migr_events,classic_merge = mass_migration(time_merge_list,pops,param1,param2)
    mass_events_list.append(mass_migr_events)
    #mass_events_list = mass_migr_events
    merge = new_arg_merge_string(pops,time_merge_list)
    
    merge_range.append(merge)
    return merge_range,mass_events_list

def merge_configuration(timestring,pops):
    merge_range = []
    merge_times = ''
    mass_events_list = []
    if event_factor != '':
        
        if 'oldest_merge' in event_factor or 'bottle_oldest_merge' in  event_factor:
            merge_times = event_factor.split('=')[1].split('-')
            steps = int(merge_times[2])
            for ti in range(int(merge_times[0]),int(merge_times[1])+steps,steps):
                num = 0
                time_merge_list = []
                merge_range,mass_events_list = process_merge(timestring,pops,merge_range,ti,time_merge_list,num,'oldest',mass_events_list)
            
        elif 'recent_merge' in event_factor or 'bottle_recent_merge' in event_factor:
            #print('merge and bottleneck')
            merge_times = event_factor.split('=')[1].split('-')
            steps = int(merge_times[2])
            for ti in range(int(merge_times[0]),int(merge_times[1])+steps,steps):
                num = 0
                time_merge_list = []
                merge_range,mass_events_list = process_merge(timestring,pops,merge_range,ti,time_merge_list,num,'recent',mass_events_list)
        else:
            num = 0
            time_merge_list = []
            merge_range,mass_events_list = process_merge(timestring,pops,merge_range,'',time_merge_list,num,'',mass_events_list)
            
    else:
        merge_range = ['']        
    return merge_range,mass_events_list



def mig_config(number_of_pops):
    migration = []
    if migration_factor != "":
        mig_value = migration_factor.split('=')
        #print(mig_value)
        mig_pops = mig_value[0]
        mig_val = mig_value[1].split('/')
       
        mig_string = '-mig ""'
        for m in mig_val:
            mig_string = mig_string +' -mig "{}={}"'.format(mig_pops,m)
        migration.append(mig_string)
        

    else:
        mig_string = ' -mig ""'
        migration.append(mig_string)
    return migration

def find_empty(param1,param2):
        if param1=='':
            return param2
        else:
            return param1


input_pars = []
bottlenecks = bottleneck(bottle_factor)
migration = mig_config(number_of_pops)
merge_t,mass_events_list= merge_configuration(timestring,number_of_pops)

simple_parameters = ''
for var in variables.keys():
    check = complicate_not_exist(var)
    
    if (not var == param1) and (not var == param2) and (check !='complicated'):
        if var not in ["recent_merge","oldest_merge","bottleneck","migration"]:

            simple_parameters = simple_parameters +' '+ arg_string(var, default[var])

if param1!='' and param2!='':
    if len(merge_t) >1:
        print('mass migr allakse apo range')
        mass_migr_events = mass_events_list
        for b in bottlenecks:
            for m in migration:
                for position,mass in enumerate(mass_migr_events):
                    middle_parameters = ''
                    for i in variables[param1]:
                        for j in variables[param2]:
                            parameters = ''
                            if complicate_par(bottlenecks,migration,merge_t):#is_merge(i,j):
                                if check_params(param1,param2):
                                    middle_parameters += b+m+mass[0]+merge_t[position] 
                                    
                                elif param1 =='recent_merge' or param1 =='oldest_merge' or param1 =='bottleneck' or param1 =='migration':
                    
                                    middle_parameters = middle_parameters + ' '+arg_string(param2, j)
                                    if j==variables[param2][-1]:
                                        middle_parameters =middle_parameters + b+m+mass[0]+merge_t[position]
                                elif param2 =='recent_merge' or param2 =='oldest_merge' or param2 =='bottleneck' or param2 =='migration':
                                    
                                    middle_parameters = middle_parameters + ' '+arg_string(param1, i)
                                    if i==variables[param1][-1]:
                                        middle_parameters = middle_parameters + b+m+mass[0]+merge_t[position]
                                
                            
                    
                    if complicate_par(bottlenecks,migration,merge_t):
                        
                        parameters = parameters+ simple_parameters+middle_parameters
                        
                        input_pars.append(parameters)

        if complicate_par(bottlenecks,migration,merge_t) == False:                 
            parameters = parameters+ simple_parameters+middle_parameters
            input_pars.append(parameters)





    else:
        #print('mass migr change')
        mass_migr_events,classic_merge = mass_migration(timestring,number_of_pops,param1,param2)

        
        p1 = -1
        p2 = -1
        seen_set1 = set([])
        seen_set2 = set([])
        
        
        parrr = ''
        for i in variables[param1]:
            if param1 not in ["recent_merge","oldest_merge","bottleneck","migration"]:

                parrr+= arg_string(param1, i)
        for j in variables[param2]:
            if param2 not in ["recent_merge","oldest_merge","bottleneck","migration"]:

                parrr+= arg_string(param2, j)
       
        


        for b in bottlenecks:
            for m in migration:
                for mass in mass_migr_events:
                    for mer in merge_t:
                        middle_parameters = ''
                        par1 =''
                        par2 = ''
                        
                        parameters = ''
                        
                        if complicate_par(bottlenecks,migration,merge_t):
                            
                            if check_params(param1,param2):
                                
                                middle_parameters = middle_parameters +b+m+mass+mer #if param1 and param2 are complicated
                            else:
                                
                                middle_parameters = parrr+ b+m+mass+mer
                                

                        else:

                            #print('got in here')
                            middle_parameters = parrr+b+m+mass+mer

                        if complicate_par(bottlenecks,migration,merge_t):
                            parameters = parameters+simple_parameters +middle_parameters+classic_merge
                            input_pars.append(parameters)

        if complicate_par(bottlenecks,migration,merge_t) == False:                 
            parameters = parameters+simple_parameters+ middle_parameters+classic_merge
            input_pars.append(parameters)
        
    
else:
    changed_parameter = find_empty(param1,param2)
    #print('nonempty',changed_parameter)

    if len(merge_t) >1:
        #print('mass migr allakse apo range')
        mass_migr_events = mass_events_list
        for b in bottlenecks:
            for m in migration:
                for position,mass in enumerate(mass_migr_events):
                    for i in variables[changed_parameter]:
                        # for j in variables[param2]:
                        parameters = ''
                        if  complicate_not_exist(changed_parameter)=='complicated':
                            
                            parameters = parameters +simple_parameters+b+m+mass[0]+merge_t[position]
                        
                        
                    input_pars.append(parameters)





    else:
        #print('mass migr paraxthike mono toy')
        mass_migr_events,classic_merge = mass_migration(timestring,number_of_pops,param1,param2)

        
        for b in bottlenecks:
            for m in migration:
                for mass in mass_migr_events:
                    for mer in merge_t:
                        
                        middle_parameters = ''
                        for i in variables[changed_parameter]:
                            parameters = ''
                            
                            if  complicate_not_exist(changed_parameter)=='complicated':
                                #print(changed_parameter,'no range sto merge')
                                parameters = parameters +b+m+mass+mer
                                
                            else:
                                middle_parameters = middle_parameters + arg_string(changed_parameter, i)
                                if i==(variables[changed_parameter][-1]):
                                    
                                    parameters = simple_parameters+ middle_parameters+mass+m+classic_merge
                                    input_pars.append(parameters)
                                    

                            check = complicate_not_exist(changed_parameter)
                            if check == 'complicated':
                                parameters = parameters+simple_parameters+classic_merge
                                input_pars.append(parameters)

#for inp in input_pars:
    #print(inp,'\n')

now = datetime.now()

output_namefile = param1+param2+'_'+str(now.strftime("%Y-%m-%d_%H:%M:%S"))
print(output_namefile,'OUT')




def sending_events(param1,param2):
    event_dictionary ={ 
        'recent_merge/bottleneck':'-event "bottle_recent_merge"',
        'bottleneck/recent_merge':'-event "bottle_recent_merge"',
        'oldest_merge/bottleneck':'-event "bottle_oldest_merge"',
        'bottleneck/oldest_merge':'-event "bottle_oldest_merge"',
        'oldest_merge':'-event "oldest_merge"',
        'recent_merge':'-event "recent_merge"',
        'bottleneck':'-event "bottleneck"'



    }
    complicated_events = ['recent_merge','oldest_merge','bottleneck']
    if (param1 in complicated_events) and (param2 in complicated_events):
        return event_dictionary[param1+'/'+param2]
    elif param1 in complicated_events:
        return event_dictionary[param1]
    elif param2 in complicated_events:
        return event_dictionary[param2]
    else:
        return '-event "noevent"'
    
event_of_run = sending_events(param1,param2)

for p in input_pars:
        os.system('python3 CLASS.py '+" "+output_namefile+" "+event_of_run+" "+p)


os.system('python3 score_plotting.py '+output_namefile)
        