import atomica as at
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import hcv_multi as hcv
import seaborn as sns
import os 
if not os.path.exists(hcv.root/"outputs"): os.makedirs(hcv.root/"outputs")
if not os.path.exists(hcv.root/"outputs"/"calibrations"): os.makedirs(hcv.root/"outputs"/"calibrations")

def calibration_plots(country:str, P, result):
    
    """ Produces a combined pdf of calibration plots
    """

    # Output name and PDF initiation
    filename = hcv.root/"outputs"/"calibrations"/f"{country}_{result.name}.pdf"
    pp = PdfPages(filename)

    # Required data and settings
    sns.set_theme(style="whitegrid")
    pops = list(P.data.pops.keys())
    pops.remove('Overall')
    plot_data = P.parsets[0]
    t_vals = np.arange(hcv.sim_start, hcv.sim_end, 1)

    # Figure 1: Population size
    fig=plt.figure(figsize=(10,10))

    for idx,pop in enumerate(pops):
        ax=fig.add_subplot(int(np.ceil(len(pops)/2)), 2, idx+1)
        ax.scatter(plot_data.get_par("alive").ts[pop].t, plot_data.get_par("alive").ts[pop].vals, s=50, c='k', label="Data")
        ax.plot(t_vals, at.PlotData(result, outputs="alive", pops=pop, t_bins=1).series[0].vals,
                  color="red", label="Model Projection")  # Model
        ax.set_title(pop, fontsize=20)
        ax.set_xlabel("Years", fontsize=20)
        ax.set_ylabel("Population", fontsize=20)
        ax.set_ylim(bottom=0)
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=18)

        # if idx==0:
        #     ax.legend(loc="best", fontsize=20)

    fig.suptitle(f"{country} Population: Model Projections", fontsize=25, style="normal", color="blue")
    fig.tight_layout()
    pp.savefig(fig)  # Attaches plot to PDF output

    # Figure 2: Prevalence
    fig=plt.figure(figsize=(10,10))

    for idx,pop in enumerate(pops):
        ax=fig.add_subplot(int(np.ceil(len(pops)/2)), 2, idx+1)
        ax.scatter(plot_data.get_par("inf").ts[pop].t, plot_data.get_par("inf").ts[pop].vals, s=50, c='k', label="Data")
        ax.plot(t_vals, at.PlotData(result, outputs="inf", pops=pop, t_bins=1).series[0].vals,
                  color="red", label="Model Projection")  # Model
        ax.set_title(pop, fontsize=20)
        ax.set_xlabel("Years", fontsize=20)
        ax.set_ylim(bottom=0)
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=18)

        # if idx==0:
        #     ax.legend(loc="best", fontsize=20)

    fig.suptitle(f"{country} Prevalence", fontsize=25, style="normal", color="blue")
    fig.tight_layout()

    # Figure 3: Tests & Treatments
    fig=plt.figure(figsize=(10,10))
    
    ax=fig.add_subplot(2, 2, 1)
    ax.scatter(plot_data.get_par("ab_notif_tot").ts["Overall"].t, plot_data.get_par("ab_notif_tot").ts["Overall"].vals, s=50, c='k', label="Data")
    ax.plot(t_vals, at.PlotData(result, pops=[{'total': pops}],
                    outputs='ab_all_m', pop_aggregation='sum', t_bins=1).series[0].vals,
              color="red", label="Model Projection")  # Model
    ax.set_title("Number of Ab+ tests", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    
    ax=fig.add_subplot(2, 2, 2)
    ax.plot(t_vals, at.PlotData(result, pops=[{'total': pops}],
                    outputs='pcr_m', pop_aggregation='sum', t_bins=1).series[0].vals,
              color="red", label="Model Projection")  # Model
    ax.set_title("Number of RNA+ tests", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    
    ax=fig.add_subplot(2, 2, 3)
    ax.scatter(plot_data.get_par("treat_tot").ts["Overall"].t, plot_data.get_par("treat_tot").ts["Overall"].vals, s=50, c='k', label="Data")
    ax.plot(t_vals, at.PlotData(result, pops=[{'total': pops}],
                    outputs='tx_m', pop_aggregation='sum', t_bins=1).series[0].vals,
              color="red", label="Model Projection")  # Model
    ax.set_title("Number of treatments", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    
    ax=fig.add_subplot(2, 2, 4)
    ax.scatter(plot_data.get_par("lt_total").ts["Overall"].t, plot_data.get_par("lt_total").ts["Overall"].vals, s=50, c='k', label="Data")
    ax.plot(t_vals, at.PlotData(result, pops=[{'total': pops}],
                    outputs='lt_m', pop_aggregation='sum', t_bins=1).series[0].vals,
              color="red", label="Model Projection")  # Model
    ax.set_title("Number of liver transplants", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)

    fig.suptitle(f"{country} Testing & Treatments", fontsize=25, style="normal", color="blue")
    fig.tight_layout()
    pp.savefig(fig)   

    # Figure 4: Incidence
    fig=plt.figure(figsize=(10,10))
    
    ax=fig.add_subplot(2, 2, 1)
    ax.stackplot(t_vals, np.vstack([at.PlotData(result, pops=pops,
                    outputs='inci_m', t_bins=1).series[i].vals for i,_ in enumerate(pops)]), labels=pops)  # Model
    ax.scatter(plot_data.get_par("hcv_inci_total").ts["Overall"].t, plot_data.get_par("hcv_inci_total").ts["Overall"].vals, s=50, c='k')
    ax.set_title("HCV incidence (new infections)", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.legend(loc='upper right')
    
    ax=fig.add_subplot(2, 2, 2)
    ax.scatter(plot_data.get_par("dc_inci_total").ts["Overall"].t, plot_data.get_par("dc_inci_total").ts["Overall"].vals, s=50, c='k', label="Data")
    ax.plot(t_vals, at.PlotData(result, pops=[{'total': pops}],
                    outputs='dc_inci', pop_aggregation='sum', t_bins=1).series[0].vals,
              color="red", label="Model Projection")  # Model
    ax.set_title("DC incidence", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    
    ax=fig.add_subplot(2, 2, 3)
    ax.scatter(plot_data.get_par("hcc_inci_total").ts["Overall"].t, plot_data.get_par("hcc_inci_total").ts["Overall"].vals, s=50, c='k', label="Data")
    ax.plot(t_vals, at.PlotData(result, pops=[{'total': pops}],
                    outputs='hcc_inci', pop_aggregation='sum', t_bins=1).series[0].vals,
              color="red", label="Model Projection")  # Model
    ax.set_title("HCC incidence", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    
    ax=fig.add_subplot(2, 2, 4)
    ax.scatter(plot_data.get_par("deaths_hcv_total").ts["Overall"].t, plot_data.get_par("deaths_hcv_total").ts["Overall"].vals, s=50, c='k', label="Data")
    ax.plot(t_vals, at.PlotData(result, pops=[{'total': pops}],
                    outputs='deaths_hcv', pop_aggregation='sum', t_bins=1).series[0].vals,
              color="red", label="Model Projection")  # Model
    ax.set_title("Number of HCV-related deaths", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)

    fig.suptitle(f"{country} Incidence & Mortality", fontsize=25, style="normal", color="blue")
    fig.tight_layout()
    pp.savefig(fig)   

    # Figure 5: Stacked outputs
    fig=plt.figure(figsize=(10,5))
    
    ax=fig.add_subplot(1, 2, 1)
    ax.stackplot(t_vals, np.vstack([at.PlotData(result, pops=pops,
                    outputs='chronic', t_bins=1).series[i].vals for i,_ in enumerate(pops)]), labels=pops)
    ax.scatter(plot_data.get_par("chronic_total").ts["Overall"].t, plot_data.get_par("chronic_total").ts["Overall"].vals, s=50, c='k')
    ax.set_title("Number of PLHCV", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    
    ax=fig.add_subplot(1, 2, 2)
    ax.stackplot(t_vals, np.vstack([at.PlotData(result, pops=pops,
                    outputs='tx_m', t_bins=1,accumulate='integrate').series[i].vals for i,_ in enumerate(pops)]), 
                 labels=pops)
    ax.set_title("Number cured", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.suptitle(f"{country} Stacked Outputs by Population Group", fontsize=25, style="normal", color="blue")
    fig.tight_layout()    
    pp.savefig(fig)  

    # Figure 6: Stacked outputs
    fig=plt.figure(figsize=(10,5))

    ax=fig.add_subplot(1, 2, 1)
    outputs = ['prop_f0f2', 'prop_f3', 'prop_f4', 'prop_dc', 'prop_hcc']
    labels = ['Proportion in F0-F2', 'Proportion in F3', 'Proportion in DC', 'Proportion in HCC']
    col = sns.color_palette("husl", 9)
    ax.stackplot(t_vals, np.vstack([at.PlotData(result, pops=[{'total': pops}],
                    outputs=outputs, t_bins=1,pop_aggregation='average').series[i].vals for i,_ in enumerate(pops)]), 
                 labels=labels, colors = col)
    ax.scatter(plot_data.get_par("prop_f0f2").ts["18–64"].t, plot_data.get_par("prop_f0f2").ts["18–64"].vals, s=50, c='k')
    ax.set_title("Number of PLHCV", fontsize=20)
    ax.set_xlabel("Years", fontsize=20)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.suptitle(f"{country} Stacked Outputs by Disease Stage", fontsize=25, style="normal", color="blue")
    fig.tight_layout()    
    pp.savefig(fig)  

    # Figure 7: Cascade
    fig=plt.figure(figsize=(10,5))

    ax=fig.add_subplot(1, 2, 1)
    
    undiag_pc, diag_pc, treat_pc = hcv.cascade_percents(result, 2017) 
    legend_labels = ['% undiagnosed', '% diagnosed but not treated', '% treated']

    bars = []
    treat = []
    undiag = []
    diag = []
    for pop in pops:
        if pop != 'Prisoners':
            bars.append(pop)
    bars.append('Overall')
    for res in bars:
        tr = treat_pc[res]
        treat.append(tr)
        dx = diag_pc[res]
        diag.append(dx)
        udx = undiag_pc[res]
        undiag.append(udx)
    treat =  np.array(treat)
    diag =  np.array(diag)   
    undiag =  np.array(undiag) 
    
    ax.bar(bars, undiag, color='r')
    ax.bar(bars, diag, bottom=undiag, color='y')
    ax.bar(bars, treat, bottom=undiag+diag, color='g')
    ax.set_title("2017", fontsize=20)
    ax.tick_params(axis='x',rotation=45, labelsize=18)

    ax=fig.add_subplot(1, 2, 2)
    
    undiag_pc, diag_pc, treat_pc = hcv.cascade_percents(result, 2030) 
    legend_labels = ['% undiagnosed', '% diagnosed but not treated', '% treated']

    bars = []
    treat = []
    undiag = []
    diag = []
    for pop in pops:
        if pop != 'Prisoners':
            bars.append(pop)
    bars.append('Overall')
    for res in bars:
        tr = treat_pc[res]
        treat.append(tr)
        dx = diag_pc[res]
        diag.append(dx)
        udx = undiag_pc[res]
        undiag.append(udx)
    treat =  np.array(treat)
    diag =  np.array(diag)   
    undiag =  np.array(undiag) 
    
    ax.bar(bars, undiag, color='r')
    ax.bar(bars, diag, bottom=undiag, color='y')
    ax.bar(bars, treat, bottom=undiag+diag, color='g')
    ax.set_title("2030", fontsize=20)
    ax.tick_params(axis='x',rotation=45, labelsize=18)
    ax.legend(legend_labels, bbox_to_anchor=(1,1),frameon=False)

    fig.suptitle(f"{country} Care Cascade", fontsize=25, style="normal", color="blue")
    fig.tight_layout()    
    pp.savefig(fig)  
    
    pp.close()
