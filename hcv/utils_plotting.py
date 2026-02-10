from hcv import utils as ut
from hcv import atomica as at
import hcv.atomica.plotting as aplt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
import seaborn as sns
import re
import geopandas as gpd
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
import matplotlib as mpl
import sciris as sc
from matplotlib.backends.backend_pdf import PdfPages
rootdir = ut.get_project_root()

mpl.rcParams["axes.edgecolor"] = "black"
mpl.rcParams["axes.linewidth"] = 1.5


def panel_plot(
    fig,
    place,
    new_fig,
    dimensions,
    col_title_text=None,
    row_title_text=None,
):
    """Panel plot configuration for a given figure.
    
    This function modifies an existing matplotlib figure by transferring an axis to a new figure and adjusting its position and labels. It also allows for the addition of optional column and row titles.
    
    Args:
        fig (matplotlib.figure.Figure): The original figure containing the axis to be modified.
        place (int): The position of the subplot in the new figure.
        new_fig (matplotlib.figure.Figure): The new figure to which the axis will be added.
        dimensions (tuple): The dimensions of the subplot grid (rows, columns).
        col_title_text (str, optional): The text for the column title. Defaults to None.
        row_title_text (str, optional): The text for the row title. Defaults to None.
    
    Returns:
        None: This function modifies the figures in place and does not return a value.
    """
    fig.set_size_inches(20, 18)
    ax = fig.axes[0]
    ax.remove()
    ax.figure = new_fig
    new_fig.add_axes(ax)
    tmp_ax = new_fig.add_subplot(*dimensions, place)
    ax.set_position(tmp_ax.get_position())
    ax.set_ylabel(" ")
    ax.set_xlabel(" ")
    tmp_ax.remove()

    plt.close(fig)

    if col_title_text is not None:
        ax.annotate(col_title_text, xy=(10, 150), xycoords="axes points", fontsize=20)
    if row_title_text is not None:
        ax.annotate(row_title_text, xy=(-80, 30), xycoords="axes points", fontsize=20,rotation=90)


def plot_outcomes(P, result, file_name, start_year=2000, end_year=2050):
    """Plot outcomes of a population model over a specified range of years.
    
    Args:
        P (Project): The project object containing population data.
        result (Result): The result object containing simulation outcomes.
        file_name (str): The title of the plot, which will be displayed as the main title.
        start_year (int, optional): The starting year for the x-axis. Defaults to 2000.
        end_year (int, optional): The ending year for the x-axis. Defaults to 2050.
    
    Returns:
        None: The function generates and saves a plot but does not return any value.
    
    Raises:
        None: This function does not raise any exceptions.
    
    Notes:
        - The function creates multiple subplots for different age groups and populations, including general populations, people who inject drugs, and prisoners.
        - It also plots various health-related metrics such as the number of people living with hepatitis C (PLHCV), deaths due to hepatitis C, treatments, new diagnoses, new infections, and new hepatocellular carcinoma (HCC) cases.
        - The plots are adjusted for aesthetics and clarity, including setting limits on the y-axis and customizing legends.
    """
    plt.rcParams["font.size"] = 10
    # Adjust plotting settings
    aplt.settings["legend_mode"] = ""
    aplt.settings["marker_edge_width"] = 1.0

    plot_data = P.data

    new_fig = plt.figure(
        figsize=(40, 18)
    )  # Make a new figure and set size (originally 22,15)
    dimensions = (4, 3)  # Lay out the subplots # (rows, columns)
    plt.subplots_adjust(hspace=0.4, wspace=0.5)  # Make some space between the plots
    plt.suptitle(file_name, fontsize=40)

    ########## Population sizes
    warnings.filterwarnings("ignore")

    d = at.PlotData(
        result, pops=["0-9_females", "0-9_males"], outputs=["alive"], project=P
    )
    fig_demo = at.plot_series(d, axis="pops", data=plot_data)
    plt.legend(["females", "_", "males"])
    plt.title("General population 0-9 years")
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    # plt.ylabel("Population size",rotation=90)
    plt.grid(False)
    panel_plot(fig_demo[0], 1, new_fig, dimensions, row_title_text="Population sizes")

    d = at.PlotData(
        result, pops=["10-17_females", "10-17_males"], outputs=["alive"], project=P
    )
    fig_demo = at.plot_series(d, axis="pops", data=plot_data, legend_mode="best")
    plt.title("General population 10-17 years")
    plt.legend(["females", "_", "males"])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 2, new_fig, dimensions)

    d = at.PlotData(
        result, pops=["18-64_females", "18-64_males"], outputs=["alive"], project=P
    )
    fig_demo = at.plot_series(d, axis="pops", data=plot_data, legend_mode=True)
    plt.title("General population 18-64 years")
    plt.legend(["females", "_", "males"])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 3, new_fig, dimensions)

    d = at.PlotData(
        result, pops=["65+_females", "65+_males"], outputs=["alive"], project=P
    )
    fig_demo = at.plot_series(d, axis="pops", data=plot_data, legend_mode=True)
    plt.title("General population 65+ years")
    plt.legend(["females", "_", "males"])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 4, new_fig, dimensions, row_title_text="Population sizes")

    d = at.PlotData(
        result, pops=["PWID_females", "PWID_males"], outputs=["alive"], project=P
    )
    fig_demo = at.plot_series(d, axis="pops", data=plot_data, legend_mode=True)
    plt.title("People who inject drugs")
    plt.legend(["females", "_", "males"])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 5, new_fig, dimensions)

    d = at.PlotData(
        result,
        pops=["Prisoners_females", "Prisoners_males"],
        outputs=["alive"],
        project=P,
    )
    fig_demo = at.plot_series(d, axis="pops", data=plot_data, legend_mode=True)
    plt.title("People in prison")
    plt.legend(["females", "_", "males"])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 6, new_fig, dimensions)

    d = at.PlotData(
        result,
        pops=[{"total": list(P.data.pops)}],
        outputs="total_hcv",
        pop_aggregation="sum",
    )
    fig_demo = at.plot_series(d, axis="results")
    plt.scatter(
        plot_data.tdve["plhcv_total"].ts["Total"].t,
        plot_data.tdve["plhcv_total"].ts["Total"].vals,
        facecolors="none",
        edgecolors="k",
    )
    plt.autoscale()
    plt.title("Number of PLHCV")
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 7, new_fig, dimensions, row_title_text="Disease burden")

    d = at.PlotData(
        result,
        pops=[{"total": list(P.data.pops)}],
        outputs="deaths_hcv",
        pop_aggregation="sum",
    )
    fig_demo = at.plot_series(d, axis="results")
    plt.scatter(
        plot_data.tdve["deaths_hcv_total"].ts["Total"].t,
        plot_data.tdve["deaths_hcv_total"].ts["Total"].vals,
        facecolors="none",
        edgecolors="k",
    )
    plt.autoscale()
    plt.title("Number of HCV deaths")
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 8, new_fig, dimensions)

    d = at.PlotData(
        result,
        pops=[{"total": list(P.data.pops)}],
        outputs="tx_m",
        pop_aggregation="sum",
    )
    fig_demo = at.plot_series(d, axis="results")
    plt.scatter(
        plot_data.tdve["treat_total"].ts["Total"].t,
        plot_data.tdve["treat_total"].ts["Total"].vals,
        facecolors="none",
        edgecolors="k",
    )
    plt.autoscale()
    plt.title("Number of treatments")
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 9, new_fig, dimensions)

    d = at.PlotData(
        result,
        pops=[{"total": list(P.data.pops)}],
        outputs="notifications_m",
        pop_aggregation="sum",
    )
    fig_demo = at.plot_series(d, axis="results")
    plt.scatter(
        plot_data.tdve["diag_nb_total"].ts["Total"].t,
        plot_data.tdve["diag_nb_total"].ts["Total"].vals,
        facecolors="none",
        edgecolors="k",
    )
    plt.autoscale()
    plt.title("Number of new diagnoses")
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 10, new_fig, dimensions, row_title_text="Disease burden")

    d = at.PlotData(
        result,
        pops=[{"total": list(P.data.pops)}],
        outputs="inci_all_m",
        pop_aggregation="sum",
        output_aggregation="sum",
    )
    fig_demo = at.plot_series(d, axis="results")
    # plt.scatter(plot_data.tdve['cirrhosis_total'].ts['Total'].t,plot_data.tdve['cirrhosis_total'].ts['Total'].vals,
    #             facecolors='none', edgecolors='C0')
    plt.autoscale()
    plt.title("Number of new infections")
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 11, new_fig, dimensions)

    d = at.PlotData(
        result,
        pops=[{"total": list(P.data.pops)}],
        outputs="hcc_inci",
        pop_aggregation="sum",
        output_aggregation="sum",
    )
    fig_demo = at.plot_series(d, axis="results")
    # plt.scatter(plot_data.tdve['dc_total'].ts['Total'].t,plot_data.tdve['dc_total'].ts['Total'].vals,
    #             facecolors='none', edgecolors='C0')
    plt.autoscale()
    plt.title("Number of new HCC cases")
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0, top=yl[1])
    plt.xlim([start_year, end_year])
    plt.grid(False)
    panel_plot(fig_demo[0], 12, new_fig, dimensions)

    # plt.show()


def plot_outcomes2(P, result, file_name, start_year=2000, end_year=2050):
    '''
    Produces a panel plot for manual calibraiton use
    :param P: atomica project
    :param result: atomica result object
    :param dirname: directory for saving the figure
    :return: panel plot
    '''
    plt.rcParams['font.size'] = 10
    # Adjust plotting settings
    aplt.settings['legend_mode'] = ''
    aplt.settings['marker_edge_width'] = 1.0
    
    plot_data = P.data

    new_fig = plt.figure(figsize=(20,18))  # Make a new figure and set size (originally 22,15)
    dimensions = (6, 5)  # Lay out the subplots # (rows, columns)
    plt.subplots_adjust(hspace=.4, wspace=0.5)  # Make some space between the plots
    plt.suptitle(file_name, fontsize =40)

    ########## Population sizes
    warnings.filterwarnings("ignore")
    
    gpop_pops_males = ['0-9_males','10-17_males','18-64_males','65+_males']
    gpop_pops_females = ['0-9_females','10-17_females','18-64_females','65+_females']

    d = at.PlotData(result, pops=['0-9_females', '0-9_males'], outputs=['alive'], project=P)
    fig_demo = at.plot_series(d, axis='pops', data=plot_data)
    plt.legend(['females', '_', 'males'])
    plt.title('General population 0-9 years')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 1, new_fig, dimensions)

    d = at.PlotData(result, pops=['10-17_females', '10-17_males'], outputs=['alive'], project=P)
    fig_demo = at.plot_series(d, axis='pops', data=plot_data, legend_mode = 'best')
    plt.title('General population 10-17 years')
    plt.legend(['females', '_', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 6, new_fig, dimensions)

    d = at.PlotData(result, pops=['18-64_females', '18-64_males'], outputs=['alive'], project=P)
    fig_demo = at.plot_series(d, axis='pops', data=plot_data, legend_mode = True)
    plt.title('General population 18-64 years')
    plt.legend(['females', '_', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 11, new_fig, dimensions)

    d = at.PlotData(result, pops=['65+_females', '65+_males'], outputs=['alive'], project=P)
    fig_demo = at.plot_series(d, axis='pops', data=plot_data, legend_mode = True)
    plt.title('General population 65+ years')
    plt.legend(['females', '_', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 16, new_fig, dimensions)

    d = at.PlotData(result, pops=['PWID_females', 'PWID_males'], outputs=['alive'], project=P)
    fig_demo = at.plot_series(d, axis='pops', data=plot_data, legend_mode = True)
    plt.title('People who inject drugs')
    plt.legend(['females', '_', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 21, new_fig, dimensions)

    d = at.PlotData(result, pops=[{'Prison_females':['Prisoners_females','Prisoners_PWID_females'], 
                                  'Prison_males':['Prisoners_males','Prisoners_PWID_males']}], outputs=['alive'], project=P)
    data_f = np.array(P.data.tdve['alive'].ts['Prisoners_females'].vals)+np.array(P.data.tdve['alive'].ts['Prisoners_PWID_females'].vals)
    data_m = np.array(P.data.tdve['alive'].ts['Prisoners_males'].vals)+np.array(P.data.tdve['alive'].ts['Prisoners_PWID_males'].vals)
    data_t = P.data.tdve['alive'].ts['Prisoners_males'].t
    fig_demo = at.plot_series(d, axis='pops', legend_mode = True)
    plt.scatter(data_t,data_f,c='b',marker='o')
    plt.scatter(data_t,data_m,c='r',marker='o')
    plt.title('People in prison')
    # plt.legend(['females', '_', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 26, new_fig, dimensions)

    d = at.PlotData(result, pops=['0-9_females', '0-9_males'], outputs=['prevalence'], project=P)
    fig_demo = at.plot_series(d, axis='pops')
    plt.title('General population 0-9 years')
    plt.legend(['females', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=1)
    # plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 2, new_fig, dimensions)

    d = at.PlotData(result, pops=['10-17_females', '10-17_males'], outputs=['prevalence'], project=P)
    fig_demo = at.plot_series(d, axis='pops')
    plt.title('General population 10-17 years')
    plt.legend(['females', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    # plt.ylim(bottom=0,top=yl[1])
    plt.ylim(bottom=0,top=1)
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 7, new_fig, dimensions)

    d = at.PlotData(result, pops=['18-64_females', '18-64_males'], outputs=['prevalence'], project=P)
    fig_demo = at.plot_series(d, axis='pops')
    plt.title('General population 18-64 years')
    plt.legend(['females', 'males'])
    plt.xlim([start_year,end_year])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.grid(False) 
    panel_plot(fig_demo[0], 12, new_fig, dimensions)

    d = at.PlotData(result, pops=['65+_females', '65+_males'], outputs=['prevalence'], project=P)
    fig_demo = at.plot_series(d, axis='pops')
    plt.title('General population 65+ years')
    plt.legend(['females', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 17, new_fig, dimensions)

    d = at.PlotData(result, pops=['PWID_females', 'PWID_males'], outputs=['prevalence'], project=P)
    fig_demo = at.plot_series(d, axis='pops', data=plot_data)
    plt.title('People who inject drugs')
    plt.legend(['females', '_', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 22, new_fig, dimensions)

    d = at.PlotData(result, pops=['Prisoners_PWID_females', 'Prisoners_PWID_males'], outputs=['prevalence'], project=P)
    fig_demo = at.plot_series(d, axis='pops', data=plot_data)
    plt.title('People in prison')
    plt.legend(['females', '_', 'males'])
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 27, new_fig, dimensions)

    d = at.PlotData(result,pops=[{'total':list(P.data.pops)}], outputs = 'total_hcv', pop_aggregation='sum')
    fig_demo = at.plot_series(d,axis='results')
    plt.scatter(plot_data.tdve['plhcv_total'].ts['Total'].t,plot_data.tdve['plhcv_total'].ts['Total'].vals,
                facecolors='none', edgecolors='k')
    plt.autoscale()
    plt.title('Number of PLHCV')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 3, new_fig, dimensions)

    d = at.PlotData(result,pops=[{'total':list(P.data.pops)}], outputs = 'deaths_hcv', pop_aggregation='sum')
    fig_demo = at.plot_series(d,axis='results')
    plt.scatter(plot_data.tdve['deaths_hcv_total'].ts['Total'].t,plot_data.tdve['deaths_hcv_total'].ts['Total'].vals,
                facecolors='none', edgecolors='k')
    plt.autoscale()
    plt.title('Number of HCV deaths')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 8, new_fig, dimensions)

    d = at.PlotData(result,pops=[{'total':list(P.data.pops)}], outputs = 'tx_m', pop_aggregation='sum')
    fig_demo = at.plot_series(d,axis='results')
    plt.scatter(plot_data.tdve['treat_total'].ts['Total'].t,plot_data.tdve['treat_total'].ts['Total'].vals,
                facecolors='none', edgecolors='k')
    plt.autoscale()
    plt.title('Number of treatments')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 13, new_fig, dimensions)

    d = at.PlotData(result,pops=[{'total':list(P.data.pops)}], outputs = 'notifications_m', pop_aggregation='sum')
    fig_demo = at.plot_series(d,axis='results')
    plt.scatter(plot_data.tdve['diag_nb_total'].ts['Total'].t,plot_data.tdve['diag_nb_total'].ts['Total'].vals,
                facecolors='none', edgecolors='k')
    plt.autoscale()
    plt.title('Number of new diagnoses')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 18, new_fig, dimensions)

    d = at.PlotData(result,pops=[{'total':list(P.data.pops)}], outputs = 'inci_all_m', pop_aggregation='sum',output_aggregation='sum')
    fig_demo = at.plot_series(d,axis='results')
    # plt.scatter(plot_data.tdve['cirrhosis_total'].ts['Total'].t,plot_data.tdve['cirrhosis_total'].ts['Total'].vals,
    #             facecolors='none', edgecolors='C0')
    plt.autoscale()
    plt.title('Number of new infections')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 23, new_fig, dimensions)

    d = at.PlotData(result,pops=[{'total':list(P.data.pops)}], outputs = 'hcc_inci', pop_aggregation='sum',output_aggregation='sum')
    fig_demo = at.plot_series(d,axis='results')
    # plt.scatter(plot_data.tdve['dc_total'].ts['Total'].t,plot_data.tdve['dc_total'].ts['Total'].vals,
    #             facecolors='none', edgecolors='C0')
    plt.autoscale()
    plt.title('Number of new HCC cases')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_demo[0], 28, new_fig, dimensions)

    #  # Stacked PLHCV
    aplt.settings['legend_mode'] = 'together'
    stacked2 = at.PlotData([result], pops=[{'Remainder population (males)': gpop_pops_males},
                                           {'Remainder population (females)': gpop_pops_females},
                                            {'People in prison (males)': ['Prisoners_males']},
                                            {'People in prison (females)': ['Prisoners_females']},
                                            {'People who inject drugs (males)': ['PWID_males']},
                                            {'People who inject drugs (females)': ['PWID_females']}
                                            ],
                           outputs='total_hcv', project=P)
    fig_prog = at.plot_series(stacked2, axis='pops', plot_type='stacked')
    plt.title('PLHCV by population group')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_prog[0], 4, new_fig, dimensions)

    stacked1 = at.PlotData([result],
                           pops=[{'total':list(P.data.pops)}],
                           outputs=['prop_f0f2', 'prop_f3', 'prop_f4', 'prop_dc', 'prop_hcc'],
                           output_aggregation='weighted', project=P)
    fig_prog = at.plot_series(stacked1, plot_type='stacked')
    plt.title('PLHCV by disease stage')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_prog[0], 9, new_fig, dimensions)

    stacked2 = at.PlotData([result],
                           pops=[{'total':list(P.data.pops)}],
                           outputs=[{'Acute':['acute_all']},
                                    {'Undiagnosed':['undiag_all']},
                                    {'Ab+ Diagnosed':['abpos_chronic']}, 
                                    {'RNA+ Diagnosed':['pcr']},
                                    {'On treatment':['treated']}],
                                    # {'Failed treatment':['treated_fail']}],
                           pop_aggregation='sum', project=P)
    fig_prog = at.plot_series(stacked2, plot_type='stacked')
    plt.scatter(plot_data.tdve['plhcv_total'].ts['Total'].t, plot_data.tdve['plhcv_total'].ts['Total'].vals,
                facecolors='k', edgecolors='k')
    plt.title('PLHCV by cascade of care')
    plt.autoscale()
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.grid(False) 
    panel_plot(fig_prog[0], 14, new_fig, dimensions)

    stacked2 = at.PlotData([result],
                           pops=[{'Remainder population (males)': gpop_pops_males},
                                 {'Remainder population (females)': gpop_pops_females},
                                 {'People in prison (males)': ['Prisoners_males']},
                                 {'People in prison (females)': ['Prisoners_females']},
                                 {'People who inject drugs (males)': ['PWID_males']},
                                 {'People who inject drugs (females)': ['PWID_females']},
                                 ],
                           outputs='inci_all_m', 
                           pop_aggregation='sum',
                           project=P)
    fig_demo = at.plot_series(stacked2,axis='pops', plot_type='stacked')
    plt.title('Number of new infections')
    yl = plt.ylim()
    plt.ylim(bottom=0,top=yl[1])
    plt.xlim([start_year,end_year])
    plt.autoscale()
    plt.grid(False) 
    panel_plot(fig_demo[0], 19, new_fig, dimensions)
    
    # cascade of the disease in end_year
    d = at.PlotData(result, project=P, pops='total', outputs='total_hcv')
    names = ['% undiagnosed', '% diagnosed but not treated', '% treated']
    
    # Plot cascade
    pops_cc = [['PWID_females'],['PWID_males'],['Prisoners_females'],['Prisoners_males'],gpop_pops_females,gpop_pops_males]
    pops_cc_names = ['PWID (females)','PWID (males)', 
                    'Prison (females)','Prison (males)',
                    'gPop (females)', 'gPop (males)']
    for year, fig_pos in zip([2020,2030],[24,29]):
        undiag_pc = dict(map(lambda key: (key, dict()), pops_cc_names))
        diag_pc = dict(map(lambda key: (key, dict()), pops_cc_names))
        treat_pc = dict(map(lambda key: (key, dict()), pops_cc_names))
        for pops, pops_name in zip(pops_cc,pops_cc_names):
            undiag_pc_pop, diag_pc_pop, treat_pc_pop = ut.cascade_percents(result, pops, year)
            undiag_pc[pops_name] = undiag_pc_pop
            diag_pc[pops_name] = diag_pc_pop
            treat_pc[pops_name] = treat_pc_pop
            
        df_cascade_pc = pd.DataFrame([undiag_pc, diag_pc, treat_pc])
        df_cascade_pc = df_cascade_pc.T
        df_cascade_pc.columns = ['% undiagnosed','% diagnosed but not treated','% treated']
        
        fig = plt.figure(figsize = (20, 16))
        fig.add_axes([0.75, 1, 0.75, 1]) 
        fig.set_size_inches(20, 16)
      
        bars = []
        treat = []
        undiag = []
        diag = []
        for pop in undiag_pc.keys():
            bars.append(pop)
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
        
        # plot bars in stack manner
        plt.bar(bars, undiag, color='r')
        plt.bar(bars, diag, bottom=undiag, color='y')
        plt.bar(bars, treat, bottom=undiag+diag, color='g')
        plt.legend(names, bbox_to_anchor=(1,1), frameon=False)
        plt.title(f'{year} care cascade')
        if year == 2020:
            plt.xticks(ticks=[],labels=[])
        else:
            plt.xticks(rotation=45, ha='right')
        plt.grid(False) 
        panel_plot(fig, fig_pos, new_fig, dimensions)
    


    df_plot = df_cascade_pc.drop('Prisoners')
    ax = df_plot.plot.bar(stacked=True, color=['r','y','g'])
    plt.xticks(rotation=45, ha='right')
    plt.legend(df_cascade_pc.columns, bbox_to_anchor=(1,1), frameon=False)
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
    ax.set_yticklabels(['0%','20%','40%','60%','80%','100%'])
    panel_plot(fig_demo[0], 24, new_fig, dimensions)
    
    plt.show()

def plot_calibration(country, cal_folder, savedir, cal_version=None):
    """Plot calibration results for a specified country.
    
    This function generates and saves a calibration plot for a given country using the specified calibration folder and directory for saving the output. It loads the calibration data, runs the simulation, and plots the outcomes.
    
    Args:
        country (str): The name of the country for which the calibration plot is to be generated.
        cal_folder (str): The folder path where calibration data is stored.
        savedir (str): The directory path where the generated plot will be saved.
        cal_version (str, optional): The version of the calibration data to use. If None, the latest version will be used.
    
    Returns:
        None: The function saves the plot as a PNG file and does not return any value.
    
    Raises:
        FileNotFoundError: If the specified calibration folder does not exist.
        ValueError: If the country is not found in the calibration data.
    """
    from hcv.parameters import iso_to_country
    print(f"Plotting country {country}")
    P = ut.project(
        country, load_calibration=True, cal_folder=cal_folder, cal_version=cal_version
    )
    result_baseline = P.run_sim(result_name="Calibration")  # run results
    # plot_outcomes(P, result_baseline, iso_to_country[country])
    # plt.savefig(savedir / f"{country}.png")
    # print(f"plot saved: {savedir}/{country}.png")
    plot_outcomes2(P, result_baseline, iso_to_country[country])
    plt.savefig(savedir / f"{country}.png")
    print(f"plot saved: {savedir}/{country}.png")
    plt.show()
    plt.close("all")


def plot_bars_impact(scens_folder, regions=None):
    """Plot bar charts to visualize the impact of different scenarios on specified regions.
    
    This function generates bar plots using seaborn to display the point estimates of various outcomes across different scenarios and regions. It calculates the averted region data from the specified folder and optionally filters the data based on the provided regions. Error bars representing the uncertainty bounds are also included in the plots.
    
    Args:
        scens_folder (str): The path to the folder containing scenario data.
        regions (list, optional): A list of regions to filter the data. If None, all regions are included.
    
    Returns:
        None: The function displays the plots directly and does not return any value.
    
    Raises:
        ValueError: If the provided folder does not contain valid scenario data.
    
    Example:
        plot_bars_impact("path/to/scenarios", regions=["AFR", "EUR"])
    """
    df_plot = ut.calc_averted_region(scens_folder)

    # Filter for specific regions
    if regions is not None:
        df_plot = df_plot[df_plot.Region.isin(regions)]

    # Plot using seaborn FacetGrid with native error bars
    sns.set(style="whitegrid")
    g = sns.catplot(
        data=df_plot,
        x="Scenario",
        y="Point Estimate",
        hue="Scenario",
        col="Outcome",
        row="Region",
        kind="bar",
        errorbar=None,
        palette="tab10",
        height=3.5,
        aspect=1.4,
        sharex=False,
        sharey=False,
        dodge=False,
        legend=False,
        edgecolor="black",
    )

    for i_row, region in enumerate(g.row_names):
        for i_col, outcome in enumerate(g.col_names):
            ax = g.axes[i_row, i_col]
            subdata = df_plot[(df_plot.Region == region) & (df_plot.Outcome == outcome)]
            for i, row in subdata.iterrows():
                xpos = list(subdata["Scenario"]).index(row["Scenario"])
                ax.errorbar(
                    x=xpos,
                    y=row["Point Estimate"],
                    yerr=[
                        [row["Point Estimate"] - row["Lower Bound"]],
                        [row["Upper Bound"] - row["Point Estimate"]],
                    ],
                    fmt="none",
                    ecolor="black",
                    capsize=3,
                    lw=1,
                )

    # Format axes and titles
    g.set_titles("{col_name}")
    g.set_axis_labels("", "")

    for i, ax in enumerate(g.axes[:, 0]):
        if g.row_names[i].upper() == "TOP10":
            ylabel = "High Burden\nCountries"
        elif g.row_names[i].upper() == "AFR":
            ylabel = "African\nRegion"
        elif g.row_names[i].upper() == "AMR":
            ylabel = "Region of the\nAmericas"
        elif g.row_names[i].upper() == "EMR":
            ylabel = "Eastern Mediterranean\nRegion"
        elif g.row_names[i].upper() == "EUR":
            ylabel = "European\nRegion"
        elif g.row_names[i].upper() == "SEAR":
            ylabel = "South East Asia\nRegion"
        elif g.row_names[i].upper() == "WPR":
            ylabel = "Western Pacific\nRegion"
        # ax.set_ylabel(g.row_names[i].upper(),fontsize=20)
        ax.set_ylabel(ylabel, fontsize=20)

    for ax in g.axes.flat:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=18)

    n_rows, n_cols = g.axes.shape

    for row_idx, row_axes in enumerate(g.axes):
        for ax in row_axes:
            ax.tick_params(axis="both", which="major", labelsize=18)
            ax.yaxis.set_major_formatter(dynamic_unit_formatter(ax))
            if row_idx < n_rows - 1:
                ax.set_xlabel("")
                ax.set_xticklabels([])
            if row_idx > 0:
                ax.set_title("")
            if row_idx == n_rows - 1:
                ax.set_xticklabels(
                    [f"Scenario {i}" for i in np.arange(1, 7, 1)]
                )  # Comment out for ppt
            else:
                ax.set_title(ax.get_title(), fontsize=20)

    g.fig.subplots_adjust(top=0.9)
    # g.savefig(fig_folder / 'epi_panel.png', dpi=300)#, bbox_inches='tight')
    # print(f'Impact plot saved: {fig_folder}/epi_panel.png')

    plt.tight_layout()
    plt.show()


def plot_bcr_map(scens_folder):
    """Plot a Benefit-Cost Ratio (BCR) map for different scenarios.
    
    This function reads a shapefile of world countries and a CSV file containing BCR values for various countries across different scenarios. It classifies the BCR values into categories, assigns colors to these categories, and then plots a map where each country is colored according to its BCR category for each scenario.
    
    Args:
        scens_folder (Path): The folder containing the 'bcr_map.csv' file with BCR values indexed by country.
    
    Returns:
        None: The function displays a plot for each scenario but does not return any values.
    """

    # Load the shapefile (update the path to where you saved it)
    world = gpd.read_file(
        str(rootdir)
        + "data/ne_110m_admin_0_countries.shp"
    )

    # Define custom colors for specific countries
    df_colors = pd.read_csv(scens_folder / "bcr_map.csv", index_col="country")

    def classify(value):
        if value < 1:
            return "BCR < 1"
        elif (value >= 1) and (value < 2):
            return "1 ≤ BCR < 2"
        elif (value >= 2) and (value < 5):
            return "2 ≤ BCR < 5"
        elif (value >= 5) and (value < 10):
            return "5 ≤ BCR < 10"
        else:
            return "BCR ≥ 10"

    color_map = {
        "BCR < 1": "darkred",
        "1 ≤ BCR < 2": "indianred",
        "2 ≤ BCR < 5": "tomato",
        "5 ≤ BCR < 10": "lightsalmon",
        "BCR ≥ 10": "mistyrose",
    }

    for scenario in df_colors.columns:
        df = df_colors
        df["Category"] = df[scenario].apply(classify)
        df["Color"] = df["Category"].map(color_map)

        # Apply colors to the GeoDataFrame
        country_colors = dict()
        for country, row in df.iterrows():
            country_colors[country] = row["Color"]

        world["color"] = world["ISO_A3_EH"].map(country_colors).fillna("lightgrey")

        # Plot the map
        fig, ax = plt.subplots(figsize=(15, 10))
        world.plot(ax=ax, color=world["color"], edgecolor="black")

        legend_elements = [
            Patch(facecolor=color_map[label], edgecolor="black", label=label)
            for label in color_map
        ]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

        ax.set_title(f"{scenario}")
        ax.set_axis_off()
        # plt.savefig(scens_folder / 'figures' / f'map_{scenario}.png', dpi=300)
        # print(f'BCR Maps saved: {scens_folder}/figures/map_{scenario}.png')
        plt.show()
        # plt.close('all')


def plot_calibration_forest(cal_folder):
    """Plot calibration data for HCV incidence and prevalent cases.
    
    This function reads calibration data from an Excel file and generates two plots: one for HCV incidence and another for HCV prevalent cases. The plots display model outputs with confidence intervals and compare them against data from the Global Hepatitis Report 2024 and The Polaris Observatory.
    
    Args:
        cal_folder (Path): The path to the folder containing the 'calibration_validation.xlsx' file.
    
    Returns:
        None: The function displays the plots directly and does not return any values.
    """
    df_plot = pd.read_excel(
        cal_folder / "calibration_validation.xlsx", sheet_name="Plot data"
    )
    columns = [
        [
            "country",
            "inci_model_pe",
            "inci_model_lb",
            "inci_model_ub",
            "inci_ghr",
            "inci_polaris",
        ],
        [
            "country",
            "plhcv_model_pe",
            "plhcv_model_lb",
            "plhcv_model_ub",
            "plhcv_ghr",
            "plhcv_polaris",
        ],
    ]
    labels = ["HCV incidence", "HCV prevalent cases"]

    for cols, label in zip(columns, labels):
        # PLHCV first
        df = df_plot[cols]
        df.columns = ["country", "model_pe", "model_lb", "model_ub", "ghr", "polaris"]
        df = df.sort_values("model_pe", ascending=True)

        # Plotting setup
        fig, ax = plt.subplots(figsize=(6, len(df) * 0.2))
        y_positions = range(len(df))

        # Plot model estimates with 95% CI
        ax.errorbar(
            df["model_pe"],
            y_positions,
            xerr=[df["model_pe"] - df["model_lb"], df["model_ub"] - df["model_pe"]],
            fmt="o",
            color="black",
            label="Model output (95% CI)",
            capsize=4,
            markersize=5,
            zorder=4,
        )

        # Plot GHR estimates
        ax.scatter(
            df["ghr"],
            y_positions,
            color="blue",
            marker="s",
            label="World Health Organization",
            zorder=3,
        )

        # Plot Polaris estimates
        ax.scatter(
            df["polaris"],
            y_positions,
            color="red",
            marker="^",
            label="© CDA Foundation",
            zorder=3,
        )

        # Set y-ticks and labels
        ax.set_yticks(y_positions)
        ax.set_yticklabels(df["country"])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))

        # Labels and legend
        ax.set_xlabel(label + " (2022)")
        # ax.set_title(label)
        ax.legend(loc="lower right", frameon=True)
        ax.grid(True, linestyle="--", alpha=0.5)

        # Tight layout
        plt.tight_layout()
        # plt.show()
        filename = cal_folder / f'{label}.eps'
        plt.savefig(filename, format='eps')
        print(f'Saved figure: {filename}.')
        

def plot_calibration_combined(cal_folder):
    """Plot calibration data for HCV incidence and prevalent cases in a 1x2 subplot layout
    with independent sorting for each metric.
    """
    df_plot = pd.read_excel(
        cal_folder / "calibration_validation.xlsx", sheet_name="Plot data"
    )
    
    columns_list = [
        ["country", "inci_model_pe", "inci_model_lb", "inci_model_ub", "inci_ghr", "inci_polaris"],
        ["country", "plhcv_model_pe", "plhcv_model_lb", "plhcv_model_ub", "plhcv_ghr", "plhcv_polaris"],
    ]
    labels = ["Incidence", "Prevalent cases"]
    titles = ["(a) Number of hepatitis C virus infections by country.",
              "(b) Number of new hepatitis C virus infections by country."]

    # Initialize 1 row, 2 columns. 
    # Increased width to 12 to provide space for y-labels on both plots.
    fig, axes = plt.subplots(1, 2, figsize=(12, len(df_plot) * 0.25))

    for i, (cols, title, label) in enumerate(zip(columns_list, titles, labels)):
        ax = axes[i]
        
        # Prepare and sort data independently for this subplot
        df = df_plot[cols].copy()
        df.columns = ["country", "model_pe", "model_lb", "model_ub", "ghr", "polaris"]
        df = df.sort_values("model_pe", ascending=True)
        
        y_positions = range(len(df))

        # Plot model estimates with 95% CI
        ax.errorbar(
            df["model_pe"],
            y_positions,
            xerr=[df["model_pe"] - df["model_lb"], df["model_ub"] - df["model_pe"]],
            fmt="o",
            color="black",
            label="Model output (95% CI)",
            capsize=4,
            markersize=5,
            zorder=4,
        )

        # Plot GHR estimates
        ax.scatter(
            df["ghr"],
            y_positions,
            color="blue",
            marker="s",
            label="World Health Organization",
            zorder=3,
        )

        # Plot Polaris estimates
        ax.scatter(
            df["polaris"],
            y_positions,
            color="red",
            marker="^",
            label="© CDA Foundation",
            zorder=3,
        )

        # Individual axis formatting
        ax.set_title(title, loc='right', fontweight='bold', pad=15)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(df["country"])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
        ax.set_xlabel(label + " (2022)")
        ax.grid(True, linestyle="--", alpha=0.5)
        
        # Add legend to each or just one depending on preference
        ax.legend(loc="lower right", frameon=True, fontsize='x-small')

    # tight_layout handles the padding between the two independent y-axis labels
    plt.tight_layout()
    
    filename = cal_folder / 'Figure_2.pdf'
    plt.savefig(filename, format='pdf', bbox_inches='tight')
    print(f'Saved subplot figure: {filename}.')
    
    filename = cal_folder / 'Figure_2.png'
    plt.savefig(filename, format='png', bbox_inches='tight')
    print(f'Saved subplot figure: {filename}.')


def dynamic_unit_formatter(ax):
    """Formats the y-axis of a given matplotlib axis object to display values in a more readable unit (e.g., billions, millions, thousands).
    
    Args:
        ax (matplotlib.axes.Axes): The matplotlib axis object to format.
    
    Returns:
        FuncFormatter: A formatter function that formats the y-axis labels based on the maximum y-value.
    """
    y_max = ax.get_ylim()[1]

    if y_max >= 1e9:
        scale = 1e9
        suffix = "B"
        precision = ".1f"
    elif y_max >= 1e6:
        scale = 1e6
        suffix = "M"
        precision = ".1f"
    elif y_max >= 1e3:
        scale = 1e3
        suffix = "K"
        precision = ".0f"
    else:
        scale = 1
        suffix = ""
        precision = ".0f"

    def formatter(x, pos=None):
        return f"{x / scale:{precision}}{suffix}"

    return FuncFormatter(formatter)


def plot_outcomes_timeseries(scens_folder, regions=["global"], scenarios="all"):
    """Plot outcomes time series for specified regions and scenarios.
    
    This function generates a time series plot of epidemiological outcomes for given regions and scenarios. It retrieves data from specified folders and an Excel file, processes the data, and creates a series of plots showing the central estimates and confidence intervals for each outcome.
    
    Args:
        scens_folder (str): The path to the folder containing scenario data.
        regions (list, optional): A list of regions to plot. Defaults to ["global"].
        scenarios (str or list, optional): Scenarios to include in the plot. If "all", all scenarios will be plotted. Defaults to "all".
    
    Returns:
        None: Displays the plot directly.
    
    Raises:
        FileNotFoundError: If the specified Excel file or folder does not exist.
        KeyError: If the specified regions or scenarios are not found in the data.
    """
    plot_data = ut.calc_outcomes_region(scens_folder, n_samples=100, regions=regions)
    regions = list(plot_data.keys())
    epi_outcomes = list(plot_data[list(plot_data.keys())[0]])
    epi_outcomes = epi_outcomes[:-1]
    data_progbook = str(rootdir) + "data/progbook_inputs.xlsx"
    df_scenarios = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name="scenarios")
    df = df_scenarios[df_scenarios.counterfactual.str.len() > 0]
    scenario_description = dict(zip(df["scenario_name"], df["description"]))
    if scenarios == "all":
        scenarios = list(
            scenario_description.keys()
        )  # list(plot_data[regions[0]][epi_outcomes[0]])
    n_rows, n_cols = 1, len(epi_outcomes)
    counterfactuals = {1: 1, 2: 2, 3: 3, 4: 3, 5: 3}

    for region in regions:
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 3), sharex=True)
        axes = axes.flatten()
        colors_order = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        for j, par in enumerate(epi_outcomes):
            ax = axes[j]
            for scen in scenarios:
                scen_name = scenario_description[scen]
                i = int(scen.split("_")[1]) + 1
                scen_label = f"Scenario {i}"
                color_nb = int(
                    re.findall(r"\d+", scen)[0]
                )  # extract scenario number as a proxy for line color

                d = plot_data[region][par][scen_name]

                # Scenario line + ribbon
                ax.plot(
                    d["years"],
                    d["central"],
                    label=f"{scen_label}",
                    lw=4,
                    color=colors_order[color_nb],
                )
                ax.fill_between(
                    d["years"],
                    d["lb"],
                    d["ub"],
                    alpha=0.1,
                    color=colors_order[color_nb],
                )

                # Counterfactual (dashed line)
                if scen == scenarios[-1]:
                    ax.plot(
                        d["years"],
                        d["cf_central"],
                        ls="--",
                        color="gray",
                        label=f"Counterfactual {counterfactuals[i]}",
                        lw=4,
                    )
            # if any(['0' in scen for scen in scenarios]):
            ax.set_title(par, fontsize=16)
            ax.grid(False)
            ax.set_xlim(2022, 2050)
            if par in ["HCV incidence", "Vaccines administered"]:
                ax.set_ylim(bottom=0)

            if j // n_cols == n_rows - 1:
                ax.set_xlabel("Year", fontsize=16)
            if j % n_cols == 0:
                ax.set_ylabel("Number of cases", fontsize=16)
            ax.yaxis.set_major_formatter(dynamic_unit_formatter(ax))
            ax.tick_params(axis="both", which="major", labelsize=14, width=1.5)
            for spine in ax.spines.values():
                spine.set_linewidth(2)
                spine.set_edgecolor("black")
    # Unique legend (bottom center)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(scenarios) + 1 if len(scenarios) < 5 else 5,
        bbox_to_anchor=(0.5, -0.2),
        fontsize=16,
        frameon=False,
    )

    # plt.suptitle(region, fontsize=14, y=1.02)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.show()


def plot_outcomes_journal(scens_folder):
    """Plot outcomes time series in a 3x3 grid with row-based scenario groups."""
    regions = ["global"]
    # Data retrieval (assuming ut and rootdir are defined in your environment)
    plot_data = ut.calc_outcomes_region(scens_folder, n_samples=100, regions=regions)
    region = "global"
    
    # Ensure we have exactly 3 outcomes for 3 columns
    epi_outcomes = list(plot_data[region].keys())
    if len(epi_outcomes) > 3:
        epi_outcomes = epi_outcomes[:3]
        
    data_progbook = rootdir / "data" / "progbook_inputs.xlsx"
    df_scenarios = pd.read_excel(pd.ExcelFile(data_progbook), sheet_name="scenarios")
    df = df_scenarios[df_scenarios.counterfactual.str.len() > 0]
    scenario_description = dict(zip(df["scenario_name"], df["description"]))
    
    # Row configuration
    scenario_groups = [
        ['scenario_0'],
        ['scenario_1'],
        ['scenario_2', 'scenario_3', 'scenario_4']
    ]
    row_titles = [
        "(a) Scenario 1.",
        "(b) Scenario 2.",
        "(c) Scenarios 3, 4, and 5."
    ]
    
    counterfactuals = {0: 1, 1: 2, 2: 3, 3: 3, 4: 3}
    colors_order = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, axes = plt.subplots(3, 3, figsize=(18, 16), sharex=False)
    
    # Storage for legend handles to ensure no duplicates
    all_handles = {}

    for row_idx, group in enumerate(scenario_groups):
        for col_idx, par in enumerate(epi_outcomes):
            ax = axes[row_idx, col_idx]
            
            for scen in group:
                scen_name = scenario_description[scen]
                scen_num = int(re.findall(r"\d+", scen)[0])
                scen_label = f"Scenario {scen_num + 1}"
                
                d = plot_data[region][par][scen_name]

                # Main Scenario Plot
                line, = ax.plot(
                    d["years"], d["central"],
                    label=scen_label,
                    lw=3.5, color=colors_order[scen_num]
                )
                ax.fill_between(
                    d["years"], d["lb"], d["ub"],
                    alpha=0.1, color=colors_order[scen_num]
                )
                
                # Collect handles for the global legend
                if scen_label not in all_handles:
                    all_handles[scen_label] = line

                # Counterfactual (plot on every subplot, but collect handle once)
                cf_idx = counterfactuals[scen_num]
                cf_label = f"Counterfactual scenarios"
                cf_line, = ax.plot(
                    d["years"], d["cf_central"],
                    ls="--", color="gray", lw=2.5,
                    # label=cf_label
                )
                if cf_label not in all_handles:
                    all_handles[cf_label] = cf_line

            # --- Formatting ---
            if col_idx == 0:
            #     ax.set_ylabel("Number of cases", fontsize=14, fontweight='bold')
            #     # Row titles on the far left
                ax.text(-0.2, 1.05, row_titles[row_idx], transform=ax.transAxes, 
                        fontsize=18, fontweight='bold', va='bottom')

            # if row_idx == 0:
            # ax.set_title(par, fontsize=18, pad=25, fontweight='bold')
            ax.set_ylabel(par, fontsize=16)
            # if row_idx == 2:
            ax.set_xlabel("Year", fontsize=16)

            ax.set_xlim(2020, 2050)
            if par in ["HCV incidence", "Vaccines administered"]:
                ax.set_ylim(bottom=0)
            
            ax.yaxis.set_major_formatter(dynamic_unit_formatter(ax))
            ax.tick_params(axis="both", labelsize=12)
            ax.grid(True, linestyle=":", alpha=0.6)
            
            for spine in ax.spines.values():
                spine.set_linewidth(1.5)

    # --- Global Shared Legend ---
    # Sort handles so Scenarios 1-5 appear in order, followed by Counterfactuals
    sorted_labels = sorted([k for k in all_handles.keys() if "Scenario" in k], 
                           key=lambda x: int(x.split()[-1]))
    cf_labels = sorted([k for k in all_handles.keys() if "Counterfactual" in k])
    
    final_labels = sorted_labels + cf_labels
    final_handles = [all_handles[lbl] for lbl in final_labels]

    fig.legend(
        final_handles, 
        final_labels,
        loc="lower center", 
        ncol=4, 
        fontsize=14, 
        frameon=True,
        bbox_to_anchor=(0.5, 0.02)
    )

    # Adjust layout to make room for bottom legend and row titles
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1, hspace=0.4, wspace=0.3)
    
    # Save as EPS
    save_path = scens_folder / "Figure_3.pdf"
    plt.savefig(save_path, format='pdf', bbox_inches='tight')
    print(f"Figure saved as {save_path}")
    
    save_path = scens_folder / "Figure_3.png"
    plt.savefig(save_path, format='png', bbox_inches='tight')
    print(f"Figure saved as {save_path}")
    plt.show()
    
def plot_calibration_panel(scens_folder, cal_folder, n_samples=100):
    """Generates a multi-page calibration panel plot for various health outcomes across multiple countries.
    
    This function reads scenario data from specified folders, processes the data to extract relevant health outcomes, and creates a PDF file containing plots for each country. Each page of the PDF displays the calibration results for a set number of countries, with each outcome represented in a separate subplot.
    
    Args:
        scens_folder (str or Path): The path to the folder containing scenario data files.
        cal_folder (str or Path): The path to the folder where the calibration panel PDF will be saved.
        n_samples (int, optional): The number of samples to process for each country. Defaults to 100.
    
    Returns:
        None: The function saves the calibration panel as a PDF file and does not return any value.
    
    Raises:
        FileNotFoundError: If the specified scenario or calibration folders do not exist.
        ValueError: If the data files do not contain the expected structure or keys.
    """
    from hcv.parameters import iso_to_country
    ## Pre-requisites for plots
    s_yr = 2010
    e_yr = 2025
    alpha = 0.05  # For 95% CrI
    sns.set_theme(style="whitegrid")

    # Define countries and outcomes
    outcomes = ["total_hcv", "inci_all_m", "tx_m", "deaths_hcv_total"]
    db_equiv = {
        "total_hcv": "plhcv_total",
        "inci_all_m": "hcv_inci_total",
        "tx_m": "treat_total",
        "deaths_hcv_total": "deaths_hcv_total",
    }
    outcome_names = ["Prevalence", "Incidence", "Treatments", "Deaths"]
    countries = ut.country_scope()

    N_COUNTRIES_TOTAL = len(countries)
    N_COLS = len(outcomes)

    # --- KEY CHANGE: Define how many countries to plot per page ---
    COUNTRIES_PER_PAGE = 10  # You can adjust this number!

    filename = cal_folder / "calibration_panel.pdf"

    # --- KEY CHANGE: Open the PdfPages object *before* the loop ---
    with PdfPages(filename) as pp:

        fig = None  # Initialize figure variable

        # Loop through all countries
        for r, country in enumerate(
            countries
        ):  # r is the *global* row index (0 to 115)

            # Calculate the row index *on the current page* (0 to 9)
            local_r = r % COUNTRIES_PER_PAGE

            # --- KEY CHANGE: Create a new figure (page) ---
            if local_r == 0:
                # If a figure already exists (from the previous page), save it
                if fig is not None:
                    fig.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2.0, w_pad=1.5)
                    pp.savefig(fig)  # Save the completed page to the PDF
                    plt.close(fig)  # Close the figure to free memory

                # Determine how many rows *this specific page* will have
                # It's COUNTRIES_PER_PAGE, unless it's the last page which might be shorter
                remaining_countries = N_COUNTRIES_TOTAL - r
                rows_on_this_page = min(COUNTRIES_PER_PAGE, remaining_countries)

                # Calculate new figure size for *this page*
                fig_width = N_COLS * 4
                fig_height = (
                    rows_on_this_page * 3.5
                )  # Base height on *actual* rows for this page
                fig = plt.figure(figsize=(fig_width, fig_height))
                print(
                    f"Creating PDF page for countries {r+1} to {r + rows_on_this_page}..."
                )

            # --- (Data loading section: No changes needed) ---
            _, db = ut.return_fw_db(country)
            country_data = dict(map(lambda key: (key, []), outcomes))
            db_data = dict(map(lambda key: (key, dict()), outcomes))

            pe = sc.load(scens_folder / country / "counter_0_central_extracted.pkl")

            tvec = None
            s_i = None
            e_i = None
            for i in range(n_samples):
                data = sc.load(scens_folder / country / f"counter_0_{i}_extracted.pkl")
                if tvec is None:
                    tvec = data["_tvec"]
                    t_list = list(tvec)
                    s_i = t_list.index(s_yr + 0.5)
                    e_i = t_list.index(e_yr + 0.5)

                for out in outcomes:
                    country_data[out].append(data[(out, "Total")][s_i : e_i + 1])
                    db_data[out] = {
                        "tvec": db.tdve[db_equiv[out]].ts["Total"].t,
                        "vals": db.tdve[db_equiv[out]].ts["Total"].vals,
                    }

            x = tvec[s_i : e_i + 1] - 0.5

            for c, (out, name) in enumerate(
                zip(outcomes, outcome_names)
            ):  # c is the column index

                sample_data = np.array(country_data[out])
                central = pe[(out, "Total")][s_i : e_i + 1]
                lb = np.percentile(sample_data, alpha / 2 * 100, axis=0)
                ub = np.percentile(sample_data, (1 - alpha / 2) * 100, axis=0)

                # --- KEY CHANGE: Calculate subplot index for the *current page* ---
                # The total number of rows is rows_on_this_page
                # The current row is local_r
                panel_idx = local_r * N_COLS + c + 1  # Matplotlib index is 1-based
                ax = fig.add_subplot(rows_on_this_page, N_COLS, panel_idx)

                # --- (Plotting section: No changes needed) ---
                ax.fill_between(x, lb, ub, alpha=0.2, color="C0", label="95% CI")
                ax.plot(x, central, lw=2, color="k", label="Central Estimate")
                ax.scatter(
                    db_data[out]["tvec"],
                    db_data[out]["vals"],
                    color="r",
                    marker="o",
                    s=20,
                    label="Data",
                )
                ax.set_xlim([s_yr, e_yr])
                ax.set_ylim([0, max(max(db_data[out]["vals"]), max(ub)) * 1.2])

                # --- KEY CHANGE: Use local_r for logic ---
                # Set titles only for the top row *of the page*
                if local_r == 0:
                    ax.set_title(f"{name}", fontsize=12)
                    if c == N_COLS - 1:
                        ax.legend(loc="best", fontsize=10)

                # Set country labels only for the left column *of the page*
                if c == 0:
                    ax.set_ylabel(
                        f"{iso_to_country[country]}",
                        fontsize=12,
                        labelpad=10,
                        rotation=90,
                        ha="right",
                    )

                ax.tick_params(axis="x", rotation=45)

                # Hide x-axis labels for all but the bottom row *of the page*
                if local_r < rows_on_this_page - 1:
                    ax.set_xticklabels([])

            # --- (Removed panel_idx increment, it's now calculated in the inner loop) ---

        # --- KEY CHANGE: Save the *final* page after the loop finishes ---
        if fig is not None:
            fig.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2.0, w_pad=1.5)
            pp.savefig(fig)
            plt.close(fig)

    # --- (Removed old single-figure save logic) ---
    print(f"Multi-page panel plot saved to: {filename}")
