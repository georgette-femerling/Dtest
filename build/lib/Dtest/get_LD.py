# Import libraries
import moments, moments.LD
import numpy as np
import demes
from collections import defaultdict
import scipy.stats as stats


def get_LD(graph,r_bins,u = 1.25e-8, normalize = True, norm_pop_idx = 0, sample_times = None, sampled_demes = None):
     '''Calculates LD per alive deme/population normalizing by itself given a demes demographic model'''
     LD_dictionary = defaultdict(list)
     if sampled_demes is None:
          alive = [deme.name for deme in graph.demes if deme.end_time == 0 ]
          sampled_demes=alive
     #print("sampling from demes:" + ",".join(sampled_demes))
     if sample_times is None:
          sample_times = np.zeros(len(sampled_demes))

     # ancestral = max(graph.demes, key=attrgetter("start_time"))
     # rhos = 4 * ancestral.epochs[0].start_size * np.array(r_bins)
     # print(ancestral.epochs[0].start_size)
     
     for i,deme in enumerate(sampled_demes):
          
          y = moments.Demes.LD(graph, sampled_demes=[deme], r=r_bins, sample_times=[sample_times[i]],u=u)
        
          # stats are computed at the bin edges - average to get midpoint estimates
          y = moments.LD.LDstats(
          [(y_l + y_r) / 2 for y_l, y_r in zip(y[:-2], y[1:-1])] + [y[-1]],
          num_pops=y.num_pops,
          pop_ids=y.pop_ids)
         
          if normalize:
               sigma = moments.LD.Inference.sigmaD2(y,normalization=norm_pop_idx)
          else:
               sigma = y
        
          LD_dictionary[deme].append(sigma)

     return LD_dictionary

def Dstat_sliced(graph,sampled_demes,r_bins,normalize = True, norm_pop_idx = 0,sample_times = None):
     if sample_times is None:
          sample_times = np.zeros(len(sampled_demes))
     #y = moments.Demes.LD(graph, sampled_demes=sampled_demes, sample_times=sample_times,rho = rhos)
     print(sampled_demes)
     rhos = 4 * graph[sampled_demes[0]].epochs[0].start_size * np.array(r_bins)
     y = moments.Demes.LD(graph, sampled_demes=sampled_demes, rho=rhos)

    # stats are computed at the bin edges - average to get midpoint estimates
     y = moments.LD.LDstats(
        [(y_l + y_r) / 2 for y_l, y_r in zip(y[:-2], y[1:-1])] + [y[-1]],
        num_pops=y.num_pops,
        pop_ids=y.pop_ids)

     if normalize:
          sigma = moments.LD.Inference.sigmaD2(y,normalization=norm_pop_idx)
     else:
          sigma = y
     return sigma

## Stats signal 
def slice_timepoints(g,time_points,yaml_filename=None):
    sliced_dict = defaultdict()
    for time in time_points:
        slicedg = moments.Demes.DemesUtil.slice(g, time)
        sliced_dict[time] = slicedg
        # option to save to yaml
        if yaml_filename != None:
            demes.dump(slicedg, yaml_filename+"_"+str(int(time))+".tmp", format='yaml', simplified=True)
    
    return sliced_dict

def get_LD_from_sliced_demes(sliced_dict,r_bins,normalize = True):
    LD_dictionary = defaultdict(list)
    for sliced in sliced_dict:
        alive = [deme.name for deme in sliced_dict[sliced].demes if deme.end_time == 0 ]
        sizes = [deme.epochs[0].start_size for deme in sliced_dict[sliced].demes if deme.end_time == 0 ]
        norm_idx = sizes.index(max(sizes))
        sigma = Dstat_sliced(sliced_dict[sliced], sampled_demes=alive, r_bins = r_bins, normalize = normalize, norm_pop_idx=norm_idx)
        for deme,i in zip(alive,range(len(alive))):
            DD = 'DD_'+str(i)+"_"+str(i)
            Dz = 'Dz_'+str(i)+"_"+str(i)+"_"+str(i)
            pi = 'pi2_'+str(i)+"_"+str(i)+"_"+str(i)+"_"+str(i)
            sigmapop = sigma.LD()[:,[sigma.names()[0].index(stat) for stat in [DD,Dz,pi]]]
            LD_dictionary[deme].append(sigmapop)
    return sigma,LD_dictionary

def calculate_signal(data_means,model,absolute=True):
    '''Calculates the Mean of the absolute log differences between the two LD decay arrays'''
    LD1=data_means
    LD2=model
    log_diff = np.log(LD1) - np.log(LD2) # calculate difference
    if absolute:
        log_diff = np.abs(log_diff)
    mld = np.mean(log_diff,axis=0)
    return mld #,np.abs(log_diff)

def calculate_signal_Dz(LDpop1,LDpop2):
    '''Calculates the Mean log difference between the two decays of D2'''
    x = LDpop1[:,1]
    y = LDpop2[:,1]
    log_diff = np.log(x) - np.log(y) # calculate difference
    mld = np.mean(np.abs(log_diff))
    return mld

def calculate_signal_D2(LDpop1,LDpop2):
    '''Calculates the Mean log difference between the two decays of D2'''
    x = LDpop1[:,0]
    y = LDpop2[:,0]
    log_diff = np.log(x) - np.log(y) # calculate difference
    mld = np.mean(np.abs(log_diff))
    return mld

def chi2_test(obs,exp,varcovs):
    '''Tests the overall significance between expected and observed values correcting by the variance of the observed data.
    Degrees of fredom are #obs - 1.
    Returns a tuple of sum(X2) and p_value'''
    chi_square = np.sum(np.square(obs - exp)/[ np.diagonal(b) for b in varcovs ],axis=0)
    p_value = 1 - stats.chi2.cdf(chi_square, (len(obs)-1))
    return(chi_square,chi_square/(len(obs)-1),p_value)

def chi2_test_perbin(obs,exp,varcovs):
    '''Tests the significance between expected and observed values correcting by the variance of the observed data.
    Returns two lists with the per-bin X2 and p_value'''
    chi_square = (np.square((obs - exp))/np.diagonal(varcovs))
    p_value = 1 - stats.chi2.cdf(chi_square, 1)
    return(chi_square,p_value)

# LD from VCF
# pops = ["pop_A", "pop_B"]
# ld_stats = moments.LD.Parsing.compute_ld_statistics(
#     vcf_path,
#     r_bins=r_bins,
#     rec_map_file=map_path,
#     pop_file=pop_file_path,
#     pops=pops
# )
