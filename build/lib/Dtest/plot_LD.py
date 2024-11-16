# Import libraries
import matplotlib.pylab as plt
import numpy as np

def plot_iicr(iirc,T): 
    plt.step(T, iirc)
    plt.xticks(fontsize= 12)
    plt.yticks(fontsize= 12)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("time ago (years)",fontsize = 14)
    plt.ylabel(r"IICR ($\it{Ne}$)",fontsize = 14)
    
def plot_comparison(LDpop1,LDpop2,labels=["Original","Size Change"],color="blue",save=False,rhos=np.logspace(-2, 2, 21)):
    
    # plot D2
    f = plt.figure(figsize=(10,5))
    ax = f.add_subplot(121)
    ax2 = f.add_subplot(122)

    ax.plot(rhos,LDpop1[:,0],label=labels[0], linewidth=3,color=color)
    ax.plot(rhos,LDpop2[:,0],'--',label=labels[1],linewidth=3,color=color)
    ax.legend()
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel(r"$\sigma_d^2$",fontsize=13)
    ax.set_xlabel(r"$\rho$",fontsize=13)

    # plot DZ
    ax2.plot(rhos,LDpop1[:,1],label=labels[0], linewidth=3,color=color)
    ax2.plot(rhos,LDpop2[:,1],'--',label=labels[1], linewidth=3,color=color)
    ax2.legend()
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylabel(r"$\sigma_{Dz}$",fontsize=13)
    ax2.set_xlabel(r"$\rho$",fontsize=13)

    plt.tight_layout()
    if save:
        plt.savefig(str(labels[0])+".Dstats.jpg",format='jpg',transparent = False)
    else:
        plt.show()
   
def plot_LD(LD_sigma, times_dic, ancestral, rhos = np.logspace(-2, 2, 21), plot_file = None, figsize = (10,20)):
    fig = plt.figure(constrained_layout=True, figsize = figsize)
    subfigs = fig.subfigures(nrows=len(LD_sigma.keys()), ncols=1)

    for pop,subfig in zip(LD_sigma,subfigs):
        subfig.suptitle(pop)
        (ax1, ax2, ax3) = subfig.subplots(nrows=1, ncols=3)
        for time_point in range(len(LD_sigma[pop])):
            ax1.plot(rhos, LD_sigma[pop][time_point][:, 0],label=str("tp_"+str(times_dic[pop][time_point])))
            ax2.plot(rhos, LD_sigma[pop][time_point][:, 1],label=str("tp_"+str(times_dic[pop][time_point])))
            ax3.plot(rhos, LD_sigma[pop][time_point][:, 2],label=str("tp_"+str(times_dic[pop][time_point])))
        
        ax1.plot(rhos, ancestral[:, 0],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax2.plot(rhos, ancestral[:, 1],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax3.plot(rhos, ancestral[:, 2],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
    
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax3.set_yscale("log")
        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")
        ax1.set_xlabel(r"$\rho$")
        ax2.set_xlabel(r"$\rho$")
        ax3.set_xlabel(r"$\rho$")
        ax1.set_ylabel(r"$\sigma_d^2$")
        ax2.set_ylabel(r"$\sigma_{Dz}$")
        ax3.set_ylabel(r"Pi2")

        ax1.legend()
        #ax2.legend()
        #ax3.legend()
    if plot_file != None :
        plt.savefig(plot_file,format='pdf',transparent = False)

def plot_LD_Linear(LD_sigma, times_dic, ancestral, rhos = np.logspace(-2, 2, 21), plot_file = None, figsize = (10,20)):
    fig = plt.figure(constrained_layout=True, figsize = figsize)
    subfigs = fig.subfigures(nrows=len(LD_sigma.keys()), ncols=1)

    for pop,subfig in zip(LD_sigma,subfigs):
        subfig.suptitle(pop)
        (ax1, ax2, ax3) = subfig.subplots(nrows=1, ncols=3)
        for time_point in range(len(LD_sigma[pop])):
            ax1.plot(rhos, LD_sigma[pop][time_point][:, 0]/ancestral[:, 0],label=str("tp_"+str(times_dic[pop][time_point])))
            ax2.plot(rhos, LD_sigma[pop][time_point][:, 1]/ancestral[:, 1],label=str("tp_"+str(times_dic[pop][time_point])))
            ax3.plot(rhos, LD_sigma[pop][time_point][:, 2]/ancestral[:, 2],label=str("tp_"+str(times_dic[pop][time_point])))
        
        ax1.plot(rhos, ancestral[:, 0],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax2.plot(rhos, ancestral[:, 1],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
        ax3.plot(rhos, ancestral[:, 2],'k--',label="Ancestral",linewidth = 1,alpha = 0.7)
    
        ax1.set_yscale("log")
        ax2.set_yscale("log")
        ax3.set_yscale("log")
        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax3.set_xscale("log")
        ax1.set_xlabel(r"$\rho$")
        ax2.set_xlabel(r"$\rho$")
        ax3.set_xlabel(r"$\rho$")
        ax1.set_ylabel(r"$\sigma_d^2/Ancestral$")
        ax2.set_ylabel(r"$\sigma_{Dz}/Ancestral$")
        ax3.set_ylabel(r"Pi2/Ancestral")

        ax1.legend()
        #ax2.legend()
        #ax3.legend()
    if plot_file != None :
        plt.savefig(plot_file,format='pdf',transparent = False)

## From moments.LD.plot - modified
def plot_ld_curves_comp(
    ld_stats,
    ms,
    vcs,
    stats_to_plot=[],
    rows=None,
    cols=None,
    statistics=None,
    fig_size=(6, 6),
    dpi=150,
    rs=None,
    numfig=1,
    cM=False,
    output=None,
    show=False,
    plot_means=True,
    plot_vcs=False,
    binned_data=True,
    ax=None,
    labels=None,
    color=None
):
    # Check that all the data has the correct dimensions
    if binned_data and (len(ld_stats.LD()) != len(rs) - 1):
        raise ValueError("binned_data True, but incorrect length for given rs.")
    if (binned_data == False) and (len(ld_stats.LD()) != len(rs)):
        raise ValueError("binned_data False, incorrect length for given rs.")


    if labels is not None:
        assert len(labels) == len(stats_to_plot)
    if labels is None:
        labels = stats_to_plot


    # set up fig and axes
    if ax is None:
        num_axes = len(stats_to_plot)
        if num_axes == 0:
            return

        if rows == None and cols == None:
            cols = len(stats_to_plot)
            rows = 1
        elif cols == None:
            cols = int(np.ceil(len(stats_to_plot) / rows))
        elif rows == None:
            rows = int(np.ceil(len(stats_to_plot) / cols))


        fig = plt.figure(numfig, figsize=fig_size, dpi=dpi)
        fig.clf()
        axes = {}


        for i, stats in enumerate(stats_to_plot):
            axes[i] = plt.subplot(rows, cols, i + 1)


    else:
        axes = [ax]


    if statistics == None:
        statistics = ld_stats.names()


    # make sure all stats are named properly
    rs = np.array(rs)
    if binned_data:
        rs_to_plot = np.array((rs[:-1] + rs[1:]) / 2)
    else:
        rs_to_plot = rs


    x_label = "$r$"
    if cM == True:
        rs_to_plot *= 100
        x_label = "cM"


    # loop through stats_to_plot, update axis, and plot
    for i, (stats, label) in enumerate(zip(stats_to_plot, labels)):
        # we don't want to plot log-scale if we predict negative values for stats
        neg_vals = False


        axes[i].set_prop_cycle(color=color)
        if plot_vcs:
            for stat in stats:
                k = statistics[0].index(stat)
                data_to_plot = np.array([ms[j][k] for j in range(len(rs_to_plot))])
                data_error = np.array(
                    [vcs[j][k][k] ** 0.5 * 1.96 for j in range(len(rs_to_plot))]
                )
                axes[i].fill_between(
                    rs_to_plot,
                    data_to_plot - data_error,
                    data_to_plot + data_error,
                    alpha=0.25,
                    label=None,
                )

        if plot_means:
            for stat in stats:
                k = statistics[0].index(stat)
                data_to_plot = np.array([ms[j][k] for j in range(len(rs_to_plot))])
                axes[i].plot(rs_to_plot, data_to_plot, "--", label=None)

        for ind, stat in enumerate(stats):
            k = statistics[0].index(stat)
            exp_to_plot = [ld_stats[j][k] for j in range(len(rs_to_plot))]
            if np.any([e < 0 for e in exp_to_plot]):
                neg_vals = True
            axes[i].plot(rs_to_plot, exp_to_plot, label=label[ind])


        axes[i].set_xscale("log")
        # don't log scale y axis for pi stats
        for stat in stats:
            if not (stat.startswith("pi2") or neg_vals):
                axes[i].set_yscale("log")


        # only place x labels at bottom of columns
        if i >= len(stats_to_plot) - cols:
            axes[i].set_xlabel(x_label)
        axes[i].legend(frameon=False, fontsize=8)
        # only place y labels on left-most column
        if i % cols == 0:
            axes[i].set_ylabel("Statistic")


    if ax is None:
        with matplotlib.rc_context(FONT_SETTINGS):
            fig.tight_layout()
            if output != None:
                plt.savefig(output)
            if show == True:
                fig.show()
            else:
                return fig