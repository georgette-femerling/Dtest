# Import libraries
import numpy as np
import demes, demesdraw
import msprime as ms
from collections import defaultdict

# Population N change model
def size_change(Ns,time_period,yaml_filename=None,plot=True,plot_filename=None,time_units="years",generation_time=29):
    
    if time_units == "generations":
        generation_time=1

    m = demes.Builder(time_units=time_units,generation_time=generation_time)
    current_time=np.array(time_period).sum()

    epochs = []
    for N,time in zip(Ns,time_period):
        epoch = dict(start_size=N,end_time=current_time-time)
        current_time=current_time-time
        epochs.append(epoch)

    m.add_deme("Deme",epochs=epochs)

    # Resolve model
    graph = m.resolve()

    # Check demographic events
    print(epochs)
    
    # option to save to yaml
    if yaml_filename != None:
        demes.dump(graph, yaml_filename, format='yaml', simplified=True)
    
    if plot:
        p = demesdraw.tubes(graph, log_time=True, num_lines_per_migration=3)
        if plot_filename != None:
            p.figure.savefig(plot_filename+".pdf")
    
    return graph 

# Cake model function
def cake_model(Ns,splits,alpha1,alpha2,time_period_merge,time_period_splits,migration_rate=1e-4,yaml_filename=None,plot=True,plot_filename=None,time_units="generations",generation_time=1):
    """
    alpha1 = proportions of previous N for split
    alpha2 = propottions of contribution to merger per splitted deme
    """
    # Check arguments
    assert len(splits) == len(alpha1) == len(time_period_splits), "Proportions and time period list must be the same length as number of split events."
    assert len(splits)+1 == len(time_period_merge), "Time period merge list must be the same length as number of split events + 1."
    if time_units == "generations":
        generation_time=1
        
    merge_events = len(splits)+1
    assert len(Ns) == merge_events, "Length of Ns list must be equal to number of split events + 1"

    total_time = np.sum(np.array(time_period_merge).sum() + np.array(time_period_splits).sum())

    m = demes.Builder(time_units=time_units,generation_time=generation_time)

    # Add first Ancestor
    m.add_deme("Ancestral",epochs=[dict(start_size=Ns[0], end_time=total_time-time_period_merge[0])])
    current_time = total_time-time_period_merge[0]
    split_b = 1
    previous = ["Ancestral"]

    #
    event = 0 
    split_i = 0
    while current_time > 0:
        if split_b:
            pops = []
            assert splits[split_i] == len(alpha1[split_i]), "Proportions list must have the same length as the number of splits"
            for pop_i,proportion in zip(np.arange(splits[split_i]),alpha1[split_i]):
                name="Split_" + str(event) + str(pop_i)
                m.add_deme(name,ancestors=previous,start_time=current_time,epochs=[dict(start_size=Ns[event]*proportion,end_time=current_time-time_period_splits[event])])
                pops.append(name)
            previous = pops
            if migration_rate > 0:
                m.add_migration(demes = pops, rate = migration_rate)
            current_time = current_time-time_period_splits[event]
            split_b = 0
            event = event + 1
        else: 
            assert len(previous) == len(alpha1[split_i]),"Length of ancestors is not equal to proportions"
            proportion = alpha2[split_i] if len(previous) > 1 else [1]
            name="Merge_" + str(event)
            m.add_deme(name,ancestors=previous,proportions=proportion,start_time=current_time,epochs=[dict(start_size=Ns[event],end_time=current_time-time_period_merge[event])])
            previous = [name]
            current_time = current_time-time_period_merge[event]
            split_b = 1
            split_i = split_i + 1

    # Resolve model
    graph = m.resolve()

    # Check demographic events
    print(graph.discrete_demographic_events()['splits'])

    # option to save to yaml
    if yaml_filename != None:
        demes.dump(graph, yaml_filename, format='yaml', simplified=True)
    
    if plot:
        p = demesdraw.tubes(graph, log_time=True, num_lines_per_migration=3)
        if plot_filename != None:
            p.figure.savefig(plot_filename+".pdf")
    
    return graph

# Load model from yaml
def load_yaml(yaml_file):
    return demes.load(yaml_file)

# https://tskit.dev/msprime/docs/stable/demography.html#sec-demography-numerical-trajectories

def _to_demes(iicr,
              T,
              Deme_name="Deme",
              yaml_filename=None,
              plot=True,time_units="years",
              generation_time=29
) -> demes.Graph: 
    """
    Takes a vector of Ns in form of iicr and T times to create a demes model
    """
    if time_units == "generations":
        generation_time=1

    Ns = np.flip(iicr)
    times = np.flip(T)

    m = demes.Builder(time_units=time_units,generation_time=generation_time)

    #current_time=np.array(time_period).sum()

    epochs = []
    for N,time in zip(Ns,times):
        epoch = dict(start_size=N,end_time=time)
        epochs.append(epoch)

    m.add_deme(Deme_name,epochs=epochs)

    # Resolve model
    graph = m.resolve()

    # option to save to yaml
    if yaml_filename != None:
        demes.dump(graph, yaml_filename, format='yaml', simplified=True)
    
    if plot:
        demesdraw.size_history(graph, log_time=True, log_size=True, title=Deme_name)
    return graph  
 
def read_relate(relate_file,generation_time=1):
    '''
    Reads relate .coal file and returns two objects: the times scaled by generation time, and the population sizes computed as 0.5/coal
    '''
    pop_size = {}
    with open(relate_file) as relatef:
        pops = relatef.readline().strip('\n').split()
        times = relatef.readline().strip('\n').split()
        times = np.array(times,dtype="float")*generation_time
        print(pops,times)
        coals = relatef.readlines()
        for line in coals:
            line = line.strip('\n').split()
            pop1 = pops[int(line[0])]
            pop2 = pops[int(line[1])]
            coal = np.array(line[2:],dtype="float")
            pop_size[(pop1,pop2)] = 0.5/coal
    relatef.close()

    return(pops,times,pop_size)

def read_psmc():
    #TODO
    return

def read_smc(smc_file,
             time_units="years",
             generation_time=29,
             yaml_filename=None
):

    assert smc_file.endswith('.csv'), 'Input must be in csv format, output from smc plot --csv'
    
    pop_sizes=defaultdict(lambda: ([], []))
    demes_dict=defaultdict()

    with open(smc_file, 'r') as f:
        next(f) #remove header
        for line in f:
            pop,time,Ne,_,_=line.strip("\n").split(',')
            pop_sizes[pop][0].append(int(float(time)))
            pop_sizes[pop][1].append(int(float(Ne)))
    
    for pop in pop_sizes:
        demes_dict[pop]=_to_demes(pop_sizes[pop][1],pop_sizes[pop][0],Deme_name=pop,time_units=time_units,generation_time=generation_time)
        # option to save to yaml
        if yaml_filename != None:
            demes.dump(demes_dict[pop], yaml_filename, format='yaml', simplified=True)
    
    return demes_dict

def get_iicr(demes_model,pop,T=None):
    """
    Returns two arrays: the Coalescence rate and the Inferred Inverse Coalescence Rate (Popsize)
    """
    m = ms.Demography.from_demes(demes_model)
    debug = m.debug()
    if np.sum(T) == None:
        T = np.concatenate([
            np.linspace(0, 1000, 2001),
            np.linspace(1000, 1.5e4, 401)[1:]
        ])
    R, _ = debug.coalescence_rate_trajectory(T, {pop: 2})
    inversed_R = 1/(2*R)

    return R,inversed_R,T

def get_iicr_croosspop(demes_model,pop1,pop2,T=None):
    """
    Returns two arrays: the Coalescence rate and the Inferred Inverse Coalescence Rate (Popsize)
    """
    m = ms.Demography.from_demes(demes_model)
    debug = m.debug()
    if T == None:
        T = np.concatenate([
            np.linspace(0, 1000, 2001),
            np.linspace(1000, 1.0e4, 401)[1:]
        ])
    R, _ = debug.coalescence_rate_trajectory(T, {pop1: 2,pop2: 2})
    inversed_R = 1/(2*R)

    return R,inversed_R,T

# Population N change model
def get_N_times_from_iicr(iicr,T):
    """
    Takes an the inverse inferred rate of coalescence and the time points and returns the times at which size changes.
    """
    if T[0] != 0:
        T.insert(0, 0)
        iicr.insert(0,iicr[0])

    assert len(T) == len(iicr), "Times and Ne have different lengths."

    previous_N = 0
    Ns = []
    times = []  
    for N,time in zip(iicr,T):
        if int(N) != int(previous_N):
            Ns.append(int(N))
            times.append(time)
            previous_N = int(N)

    return np.flip(Ns),np.flip(times)

def size_change_from_iicr(iicr,T,Deme_name="Deme",yaml_filename=None,plot=True,plot_filename=None,time_units="years",generation_time=29):
    """
    Takes a vector of Ns in form of iicr and T times to create a model
    """
    if time_units == "generations":
        generation_time=1

    Ns,times = get_N_times_from_iicr(iicr,T)

    
    m = demes.Builder(time_units=time_units,generation_time=generation_time)

    #current_time=np.array(time_period).sum()

    epochs = []
    for N,time in zip(Ns,times):
        epoch = dict(start_size=N,end_time=time)
        epochs.append(epoch)

    m.add_deme(Deme_name,epochs=epochs)

    # Resolve model
    graph = m.resolve()

    # option to save to yaml
    if yaml_filename != None:
        demes.dump(graph, yaml_filename, format='yaml', simplified=True)
    
    if plot:
        p = demesdraw.tubes(graph, log_time=True, num_lines_per_migration=3)
        if plot_filename != None:
            p.figure.savefig(plot_filename+".pdf")
    
    return graph  
    