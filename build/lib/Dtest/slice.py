# Import libraries
import moments, moments.LD
from collections import defaultdict
import numpy as np
import demes

def slice_timepoints(g,time_points,yaml_filename=None):
    sliced_dict = defaultdict()
    for time in time_points:
        slicedg = moments.Demes.DemesUtil.slice(g, time)
        sliced_dict[time] = slicedg
        # option to save to yaml
        if yaml_filename != None:
            demes.dump(slicedg, yaml_filename+"_"+str(int(time))+".tmp", format='yaml', simplified=True)
    
    return sliced_dict

def sample_timepoints(graph,anc_end_time=1000000):
    times = []
    times_dic = defaultdict()
    for deme in graph.demes:
        time = []
        for epoch in deme.epochs:
            start = (epoch.start_time if epoch.start_time != float('inf') else anc_end_time)
            midway = (start + epoch.end_time)/2
            time.extend([start-1, midway, epoch.end_time])
        times = np.concatenate([times,time])
        times_dic[deme.name] = time
    return np.unique(np.flip(times)),times_dic