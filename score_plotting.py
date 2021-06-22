import os
os.environ["OMP_NUM_THREADS"] = "10"
import numpy as np


import matplotlib
import matplotlib.pyplot as plt
import argparse
import pathlib

from matplotlib.ticker import ScalarFormatter

# Simple helper functions
def cast(t,value):
    ## Checking type of a value
    
    if t == 'int':
        return int(value)
    elif t == 'float':
        return float(value)

def get_type(header,i):
    ## Categorizes variables depending on the type of their values

    #integers=['sample','recent_merge','oldest_merge','bottleneck_pops']
    floats=['mut_rate','recomb_rate','length','migration']
    ret = 'int'
    if header[i] in floats:
        ret = 'float'
    return ret

def prettify(i):
    ## Changes name of variables that will be shown in plots

    d = {
        "migration": "Migration",
        "recomb_rate": "Recombination Rate",
        "mut_rate": "Mutation Rate",
        'length': 'Length',
        'sample': 'Number of Samples',
        'oldest_merge': 'Time of Merge',
        'recent_merge': 'Time of Merge',
        'bottleneck_pops': 'Bottleneck in pops',
        "bottle_end": 'End of Bottleneck',
        'bottle_st': 'End of Bottleneck'
    }
    return d[i]

def get_field(line,header,search_var):
    ## Gets requested field in line of a file

    position = None
    for i,m in enumerate(header):
        if m == search_var:
            position= i
            break
    return line[position]
# Returns number of common elements between two sets
def identification(a,b): 
    ## Checking if valid parameters were given
    return len(set(a) & set(b))

# Call appropriate plotting functions

def config(header):
    ## Finds variables in header of a file

    parameters = ['length','mut_rate','recomb_rate','sample']
    complicated_param = ['oldest_merge','recent_merge','bottleneck_pops','migration','pop_bottled']
    all_parameters = ['length','mut_rate','recomb_rate','sample','oldest_merge','recent_merge','bottleneck_pops','migration','pop_bottled']
    
    if identification(all_parameters,header)==2:
        return '3D'
       
    elif identification(parameters,header)==1:
        return 'simple'

    elif identification(complicated_param,header)==1:
        return 'complicated'

def get_type_of_score(name_of_file):
    ## Gets type of score

    d = {
        'pca_eucl': 'pca_euclideian',
        'fst': 'Fst',
        'pca_manh': 'pca_manhattan',
        'pca_bray':'pca_braycurtis',
        'tsne_bray':'tsne_braycurtis',
        'tsne_eucl':'tsne_euclideian',
        'tsne_manh':'tsne_manhattan'
    }
    for key in d.keys():
        if key in name_of_file:
            return d[key]

def create_plot(filename):
    ## Distinquises type of file that is about to be plotted and type of plot that will be used

    
    with open(filename,'r') as f:
        
        # if 'pca_eucl' in filename:
        #     type_of_score = 'pca_euclideian'
        # elif 'pca_manh' in filename:
        #     type_of_score = 'pca_manhattan'
        # elif 'pca_bray' in filename:
        #     type_of_score = 'pca_braycurtis'
        # elif 'fst' in filename:
        #     type_of_score = 'Fst'
        # elif 'tsne_bray' in filename:
        #     type_of_score = 'tsne_braycurtis'
        # elif 'tsne_eucl' in filename:
        #     type_of_score = 'tsne_euclideian'
        # elif 'tsne_manh' in filename:
        #     type_of_score = 'tsne_manhattan'
        type_of_score = get_type_of_score(filename)
        header = f.readline().strip().split(' ')
        param_type = config(header)
        
        if param_type == 'simple':
            
            plotting_2D(f,header,type_of_score)
        elif param_type == 'complicated':
            
            complicated_plotting_2D(f,header,type_of_score)
        elif param_type == '3D':
            
            plotting_3D(f,header,type_of_score)

# 2D plotting

def final_plotting(header,path, type_of_score, xlabel, ylabel, title, var1, var2, start_of_bottle):
    ## 2D Helper function

    i1=0
    
    if title is None:
        title = ''
    figure, ax = plt.subplots(1)
    ax.scatter(var1,var2,color='r',zorder=1,s=100)
    ax.plot(var1,var2,color='b',linewidth=4.0,zorder=2)
    if xlabel == 'Migration Rate':
        plt.xlim(0,max(var1))

    else:
        plt.xlim(0)
        
    ax.set_ylim(0,105)
    
    
    final_title = ' Verified Scenario changing {}\n'.format(prettify(header[i1]))
    final_title = final_title + title
    plt.title(final_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    plt.savefig(path)
    
def plotting_2D(f, header, type_of_score):
    ## Plot of single simple variable(e.g samples)

    
    path=pathlib.Path(os.path.realpath(f.name)).parent.parent
    path = str(path) + '/'+type_of_score+'_result_plot.png'
    var1 = []
    var2 = []
    is_float = identification(['mut_rate','recomb_rate','length'],header)
    ylabel= type_of_score+' Success Score'

    for line in f:
        line = line.split(' ')
        
       
        if is_float:
            var1.append(float(line[0]))
        else:
            var1.append(int(line[0]))
        var2.append(int(line[1]))
    
    if 'mut_rate' in header:
        xlabel = 'Mutation Rate'
    elif 'recomb_rate' in header:
        xlabel = 'Recombination Rate'
    elif 'length' in header:
        #var1 = [1e-7*i for i in var1]
        xlabel = 'Length'
    elif 'sample' in header:
        xlabel = 'Number of Samples'

    final_plotting(header,path,type_of_score,xlabel,ylabel,None,var1,var2,start_of_bottle=None)

def complicated_plotting_2D(f, header, type_of_score):
    ## Plot of single complicated variable(e.g bottleneck)
    path=pathlib.Path(os.path.realpath(f.name)).parent.parent
    path = str(path) + '/'+type_of_score+'_result_plot.png'
    
    
    var1 = []
    var2 = []
    bottleneck_val = []
    
    if ('bottleneck_pops' not in header) and ('pop_bottled' not in header):
        pos=(0,2)
    else:
        pos=(1,3)
        
    
    for line in f:
        line = line.split(' ')

        if 'migration' in header:
            var1.append(float(line[pos[0]]))
            var2.append(int(line[pos[1]]))
        elif 'oldest_merge' in header or 'recent_merge' in header:
            var1.append(int(line[pos[0]]))
            var2.append(int(line[pos[1]]))

        elif 'bottleneck_pops' in header or 'pop_bottled' in header:
            bottleneck_val.append((int(line[pos[0]]),int(line[pos[1]])))
        
        
    if len(bottleneck_val)>1:
        bottleneck_val = sorted(bottleneck_val, key=lambda tup: tup[0])
        for thesi,_ in enumerate(bottleneck_val):
            var1.append(bottleneck_val[thesi][0])
            var2.append(bottleneck_val[thesi][1])
    
    xlabel=header[pos[0]]
    ylabel=type_of_score+' Success Score'
    start_of_bottle = None
    checked_pops = []
    if 'migration' in header:

        migrated_pops = line[1].split('|')
        for th,i in enumerate(migrated_pops):
            pair = migrated_pops[th]# .split(',')
            if pair not in checked_pops:
                checked_pops.append(pair)
        title = 'migration between pop{}-pop{}'.format(migrated_pops[0][0],migrated_pops[0][2])
        
        for i in checked_pops[1:]:
            elements = i.split(',')
            title = title + ' and pop{}-pop{}'.format(elements[0],elements[1])
        xlabel = 'Migration Rate'

    elif 'oldest_merge' in header or 'recent_merge' in header:

        merged_pops = line[1]
        
        pair = merged_pops.split(',')
        title = 'There was a merge between pop{}-pop{}'.format(pair[0],pair[1])
        xlabel = 'Time of Merge'

        
    elif 'bottleneck_pops' in header or 'pop_bottled' in header:
        
        bottled_populations = line[2].split(',')
        title = 'There was bottleneck in '
        title = title + 'pop{}'.format(int(bottled_populations[0]))
        for i in bottled_populations[1:]:
            
            title = title + ' and pop{}'.format(i)

        start_of_bottle = line[0]
        title = title + '\n which started {} years ago'.format(start_of_bottle)
        
        xlabel ='End of Bottleneck'

    final_plotting(header,path,type_of_score,xlabel,ylabel,title,var1,var2,start_of_bottle)

# Heatmap-related functions

def plotting_3D(f, header, type_of_score):
    ## Heatmap helper
    
    path=pathlib.Path(os.path.realpath(f.name)).parent.parent
    path = str(path) + '/'+type_of_score+'_result_plot.png'
    
    var1 = []
    var2 = []
    var3 = []
    i1=0
    i2=1
    start_of_bottle = None
    #is_float = identification(['mut_rate','recomb_rate','length','migration'],header)!=0
    (type_1,type_2)=(get_type(header,i1), get_type(header,i2))
    
    ylabel= type_of_score+' Success Score'
    
    for line in f:
        
        line = line.split(' ')
        
        var1.append(cast(type_1,line[i1]))
        var2.append(cast(type_2,line[i2]))
        var3.append(int(get_field(line,header,'succ_score')))

    
    xlabel = prettify(header[i2])
    ylabel = prettify(header[i1])
    title = ' Verified Scenario changing {},{}\n'.format(prettify(header[i1]),prettify(header[i2]))
    
    
    if 'migration' in header:
        
        migrants = get_field(line,header,'migr_pops')
        if migrants != '0':
            migrated_pops = migrants.split('|')
            migrated_pops = np.unique(migrated_pops)
            pair = migrated_pops[0].split(',')
            title = title+ 'Migration: pop{}-pop{}'.format(pair[0],pair[1])
            for i in migrated_pops[1:]:
                elements = i.split(',')
                title = title + ' ,pop{}-pop{}'.format(elements[0],elements[1])
            

    if ('oldest_merge' in header) or ('recent_merge' in header):

        merged_pops = get_field(line,header,'merged_pops')
        
        pair = merged_pops.split(',')
        title = title+ ' \nMerge: pop{}-pop{}'.format(pair[0],pair[1])
        
        
    if ('bottleneck_pops' in header) or ('pop_bottled' in header):
        if 'bottleneck_pops' in header:
            bottled_populations = get_field(line,header,'bottleneck_pops').split(',')
        elif 'pop_bottled' in header:
            bottled_populations = get_field(line,header,'pop_bottled').split(',')


        #bottled_populations = get_field(line,header,'bottleneck_pops').split(',')
        title = title+' \nBottleneck: '
        title = title + 'pop{}'.format(int(bottled_populations[0]))
        for i in bottled_populations[1:]:
            
            title = title + ' ,pop{}'.format(i)

        start_of_bottle = get_field(line,header,'bottle_st')
        title = title + '\n Start: {} years ago'.format(start_of_bottle)
        #xlabel ='End of Bottleneck'

    
    
    fig, ax = plt.subplots()
    
    im, cbar,scores_array,x,y = heatmap(xlabel,ylabel,var3,var1,var2,vmin=0,vmax=100,cmap="YlGnBu", cbarlabel=type_of_score+' Success Score',ax=ax)
    
    texts = annotate_heatmap(im,scores_array,x,y,valfmt="{x:.1f} %",threshold=50)
    ax.set_title(title)
    fig.tight_layout()
    # ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    # ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.show()
    plt.savefig(path)


def data_to_plot(var1,var2,var3):
    ## Prepare data for heatmap
    
    
    
    diction = {}
    for i,_ in enumerate(var1):
            diction[(var1[i],var2[i])] = var3[i]
   
    
    dic = sorted(diction, key = lambda x: (x[0]),reverse=True)
    values = []
    for key in dic:
        values.append(diction[key])
    
    dimension_x = len(set(var1))
    dimension_y = len(set(var2))
    data = np.zeros((dimension_x,dimension_y))
    rows,col = data.shape
    count = 0
    
    for i in range(0,rows):
        for j in range(0,col):
            data[i,j] = int(values[count])
            count+=1

    
    
    return data

def annotate_heatmap(im,scores_array,x,y,valfmt,textcolors=["black", "white"],threshold=None, **textkw):
    ## Customize heatmap visualization
    if not isinstance(scores_array, (list, np.ndarray)):
        scores_array = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(scores_array.max())

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(scores_array.shape[0]):
        for j in range(scores_array.shape[1]):
            kw.update(color=textcolors[int(im.norm(scores_array[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(scores_array[i, j], None), **kw) 
            texts.append(text)

    return texts

def heatmap(xlabel,ylabel,var3,var1,var2,ax=None,cbar_kw={}, cbarlabel="",**kwargs):
    ## Making heatmap visualization for 3 variables 

    y = sorted(var1)
    x = sorted(var2)

    
    x = sorted(set(x))
    y= sorted(set(y),reverse=True) #setting axis
    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    scores_array = data_to_plot(var1,var2,var3)
    im = ax.imshow(scores_array, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    ax.set_xticks(np.arange(scores_array.shape[1]))
    ax.set_yticks(np.arange(scores_array.shape[0]))
    
    ax.set_xticklabels(x)
    ax.set_yticklabels(y)


    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for _, spine in ax.spines.items():
        spine.set_visible(False)

    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar,scores_array,x,y

def results_visualization(folder):
    ## Calls functions of visualization for every file in folder
    folder = folder+'/scores_to_plot/'
    count = 0

    for filename in os.listdir(folder):
        create_plot(folder+str(filename))
        
        count+=1


parser = argparse.ArgumentParser()
parser.add_argument('folder', help='folder where output needs to be plotted')  
args = parser.parse_args()
folder = args.folder
results_visualization(folder)
os.system("command -v osascript > /dev/null && osascript -e 'display notification \"PLOT Done\"'")



