#!/bin/python

"""****************************************************************************
This interactive python routine facilitates execution of DayCent simulations
and results analysis for bioenergy ecosystem dynamics assessments. It is
currently organized around the generation of Figures 3a-c for the manuscript
"Greenhouse Gas Emission Reduction from Biofuels and Ecosystem Carbon Storage",
and structured such that:
   * model execution and results analysis are de-coupled for convenience
   * ensemble runs across multiple sites are accommodated
   * analysis results are automatically archived with full provenance
To add additional sites or scenarios to the model run ensemble, create a new
site directory with weather.wth and site.100 files and add the associated
scenario details to the run_combos data structure.
****************************************************************************
"""

import os
import csv
import glob
import shutil
import subprocess
import time
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import sqlite3


# A few constants, since the normal constants.py is not necessarily accessible
n_to_n2o = 44.013/28.014
n2o_gwp100 = 298.0


def stacked_line_plot(axis_object, xs, ys_list, color_hatch_label_list, legend=False):
    previous_cum_y = []
    current_cum_y = []
    for i, ys in enumerate(ys_list):
        if i == 0:
            current_cum_y = ys
            color, hatch, label = color_hatch_label_list[i]
            axis_object.plot(xs, ys, color=color)
            previous_cum_y = ys
        else:
            current_cum_y += ys
            color, hatch, label = color_hatch_label_list[i-1]
            axis_object.plot(xs, current_cum_y, color=color)
            if hatch:
                axis_object.fill_between(xs, current_cum_y, ys_list[0], facecolor='none', hatch=hatch, edgecolor=color,
                                         zorder=-i)
                axis_object.bar([0], [0], facecolor='none', hatch=hatch, edgecolor=color, label=label)
            else:
                axis_object.fill_between(xs, current_cum_y, ys_list[0], facecolor=color, linewidth=0.0, zorder=-i)
                axis_object.bar([0], [0], color=color, label=label)
            previous_cum_y = current_cum_y
    axis_object.plot(xs, current_cum_y, color='k', label='Cumulative NPP')
    if legend:
        axis_object.legend(loc=2, frameon=False, prop={'size': 6.5})


def gm2_to_Mgha(data_list):
    for i, datum in enumerate(data_list):
        data_list[i] = float(datum) * 0.01
    return data_list


def DDcentEVI(ddc_fpath, sch_file, target_path, id, out_files, ddclist_fpath="", site_file_out=False):
    """Runs the specified versions of DDcentEVI and DDClist100 for the specified schedule
    file located in the target directory, saving (renaming) all specified daily output
    files, with all output files named as per the specified ID. Note that any .out files
    to preserve must also be specified in the outfiles.in file in referenced library
    archive.
    Args:
        ddc_fpath- file/path to UNIX-executable DailyDayCent program (str)
        ddclist_fpath- file/path to UNIX-executable list100 program (str)
        sch_file- filename of schedule file in target directory (str)
        target_path- directory where all model inputs are located and results are saved (str)
        id- name to be given to all model output files (str)
        out_files- list of daily output .out files to save (list of str)
    """
    os.chdir(target_path)
    if site_file_out:
        subprocess.call("%s -s %s -n %s -W %s.100" % (ddc_fpath, sch_file, id, id), shell=True)
    else:
        subprocess.call("%s -s %s -n %s" % (ddc_fpath, sch_file, id), shell=True)
    if ddclist_fpath:
        subprocess.call("%s %s %s outvars.txt" % (ddclist_fpath, id, id), shell=True)
    for file in out_files:
        os.rename(target_path+file, target_path+id+"_"+file)


def DDC_detail(work_path, results_path, detail_id, detail_start_year, detail_end_year):
    # r_fpath = '/data/paustian/AFRI/simulations/visualization/growthlim_diagnostics_forUnix_112414.R'
    r_fpath = '/data/paustian/GEVO_jeff_kent/autorun_testing/Rdir/growthlim_diagnostics_forUnix_022715.R'
    command = 'Rscript %s %s "" TRUE TRUE FALSE FALSE 2.5 FALSE TRUE %s detail_%s png 12 FALSE %i %i TRUE TRUE' % \
              (r_fpath, work_path, results_path, detail_id, detail_start_year, detail_end_year)
    print command
    # raw_input("Press RETURN to continue...")
    subprocess.call("%s" % command, shell=True)
    # raw_input("Press RETURN to continue...")


def auto_type(datum):
    """Convert data types to integer or float if possible
    """
    try:
        return int(datum)
    except:
        try:
            return float(datum)
        except:
            return str(datum)


def summarize_out(out_fpath, head_skip, tail_skip, avg_cols, diff_cols):
    """Reads space-delimited .lis or .out output files (including skipping headers and
    dummy head or tail data rows), and returns averages of entire columns or averages of
    annual differences within a column.
    Args:
        out_fpath- full path and filename of output file to be analyzed (str)
        head_skip- number of rows to skip at beginning of file, incl. headers and dummy data (int)
        tail_skip- numbers of rows to skip at end of file, typically dummy data (int)
        avg_cols- columns for which to return averages, counting from 0 (list of int)
        diff_cols- columns for which to return (list of int)
    """
    #skip the proper number of headers, read each subsequent line into a 2D list, assuming
    #space-delimited, and delete the 2D list initialization row and any dummy last lines
    cont = [[]]
    input = open(out_fpath, 'r')
    for i in range(head_skip):
        next(input)
    for line in input:
        cont.append(line.split())
    del cont[0]
    for j in range(tail_skip):
        del cont[-1]

    #convert all entries in the 2D list to float format where possible, or zero in the case
    #of very small numbers recorded in scientific notation
    for k in range(len(cont)):
        for l in range(len(cont[k])):
            cont[k][l] = float(cont[k][l])

    #loop through and compute each avg_cols entry
    avg_cols_out = []
    for m in avg_cols:
        tot = 0
        count = 0
        for p in range(len(cont)):
            tot += cont[p][m]
            count += 1
        avg_cols_out.append(tot/count)

    #loop through and compute each diff_cols entry
    diff_cols_out = []
    for n in diff_cols:
        tot = 0
        count = 0
        for p in range(len(cont)-1):
            tot += cont[p+1][n]-cont[p][n]
            count += 1
        diff_cols_out.append(tot/count)

    return avg_cols_out, diff_cols_out


def read_full_out(out_fpath, head_skip, tail_skip, delim='s'):
    """Reads space-delimited .lis or .out output files (including skipping headers and
    dummy head or tail data rows), and returns averages of entire columns or averages of
    annual differences within a column.
    Args:
        out_fpath- full path and filename of output file to be analyzed (str)
        head_skip- number of rows to skip at beginning of file, incl. headers and dummy data (int)
        tail_skip- numbers of rows to skip at end of file, typically dummy data (int)
    """
    # read all lines into a 2D list, delete the list initialization and specified head/tail lines
    input = open(out_fpath, 'rU')
    mylist = [[]]
    if delim == 's':
        for i in range(head_skip):
            next(input)
        for line in input:
            mylist.append(line.split())
        del mylist[0]
        for j in range(tail_skip):
            del mylist[-1]

    else:
        if delim == 't':
            lines = csv.reader(input, delimiter="\t")
        elif delim == 'c':
            lines = csv.reader(input)
        for line in lines:
            mylist.append(line)
        del mylist[0]
        for i in range(head_skip):
            del mylist[0]
        for j in range(tail_skip):
            del mylist[-1]

    # convert all entries in the 2D list to float format where possible, or zero in the case
    # of very small numbers recorded in scientific notation
    for k in range(len(mylist)):
        for l in range(len(mylist[k])):
            mylist[k][l] = auto_type(mylist[k][l])

    #return transpose
    mylist = zip(*mylist)
    return mylist


def yie_ghg_sum(id, work_path):
    """For the specified bioenergy or BAU scenario, this function takes basic output file
    column-averaged output from read_out(), and performs the appropriate unit conversions
    to report scenario biomass yield, grain yield, change in SOM, N2O emissions, and total
    CO2 equivalent GHGs.
    Args:
        id- string describing simulating, incl. soil, LUC scenario, bioenergy/BAU (str)
    """
    #call read_out() to analyze the simulation .lis and year_summary.out output
    avg_cols_out1, diff_cols_out1 = summarize_out(work_path+id+'.lis', 4, 1, [6,8], [7,9])
    avg_cols_out2, diff_cols_out2 = summarize_out(work_path+id+'_year_summary.out', 1, 0, [1], [])
    cmrvst = avg_cols_out1[0]   #gC/m2/y
    tcrem = diff_cols_out1[0]    #gC/m2/y
    cgrain = avg_cols_out1[1]   #gC/m2/y
    d_somsc = diff_cols_out1[1]   #gC/m2/y
    n2oflux = avg_cols_out2[0]   #gN2O-N/m2/y

    #unit conversions & data recording
    hbmyield = (cmrvst/0.45)/100.0  #Mg-BM/ha/y
    wbmyield = (tcrem/0.45)/100.0  #Mg-BM/ha/y
    gryield = (cgrain/0.45)/100.0  #Mg-grain/ha/y
    d_som = d_somsc/100.0  #Mg/ha/y
    n2o = n2oflux*(44.0/28.0)*10000.0*.001
    ghg = (-1*d_som)+(298.0*n2o*0.001)

    case = id.split('_')[-1]
    sum = [case, round(hbmyield,2), round(wbmyield,2), round(gryield,2), round(d_som,4),
           round(n2o,3), round(ghg,3)]
    return sum


def bioenergy_ghg_avoid_gCO2eq_m2(biomass_gC_m2_array, conversion_tech):
    """Based on 'JohnF_BiomassBioenergy_7_28_17-JLF-v2.xlsx', comparison tab
    """

    if conversion_tech not in ('current', 'mature', 'CCS'):
        print "*** ERROR: unknown conversion technology specified"
        print
        return

    # GHG intensity factors, all in units of g CO2eq (kg dry delivered biomass)-1
    ghg_int_biofuel_supplychain_g_kg = 173.20

    gross_GHG_int_avoided_g_kg = -9999999
    CCS_marg_GHG_int_g_kg = -9999999

    if conversion_tech == 'current':
        gross_GHG_int_avoided_g_kg = 705.02
        CCS_marg_GHG_int_g_kg = 0

    elif conversion_tech == 'mature':
        gross_GHG_int_avoided_g_kg = 1110.0
        CCS_marg_GHG_int_g_kg = 0

    elif conversion_tech == 'CCS':
        gross_GHG_int_avoided_g_kg = 1110.0
        CCS_marg_GHG_int_g_kg = 766.74

    # calculate individual GHG intensity terms from biomass supply array
    biomass_c_conc = 0.45
    biomass_array_kg_m2 = (biomass_gC_m2_array / biomass_c_conc) * 0.001   # kg biomass m-2 y-1

    gross_bioenergy_ghg_mitigation_array_g_m2 = gross_GHG_int_avoided_g_kg * biomass_array_kg_m2   # g CO2eq m-2 y-1
    supply_chain_ghg_array_g_m2 = ghg_int_biofuel_supplychain_g_kg * biomass_array_kg_m2   # g CO2eq m-2 y-1
    CCS_marg_ghg_array_g_m2 = CCS_marg_GHG_int_g_kg * biomass_array_kg_m2   # g CO2eq m-2 y-1
    net_bioenergy_ghg_mitigation_array_g_m2 = gross_bioenergy_ghg_mitigation_array_g_m2 + CCS_marg_ghg_array_g_m2 - supply_chain_ghg_array_g_m2

    return net_bioenergy_ghg_mitigation_array_g_m2, supply_chain_ghg_array_g_m2, CCS_marg_ghg_array_g_m2


def copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, library_path, bin_fpath=""):
    """Clear out workspace, then copy in schedule, weather, soil & library files
    """
    for file in glob.glob(os.path.join(work_path, '*.*')):
        os.remove(file)
    shutil.copy(site_fpath, work_path)
    shutil.copy(schedule_fpath, work_path)
    shutil.copy(weather_fpath, work_path)
    if bin_fpath:
        shutil.copy(bin_fpath, work_path)
    for file in glob.glob(os.path.join(library_path, '*')):
        shutil.copy(file, work_path)
    return


def clean_out(work_path, move_path, save_types):
    """Archive key results and clean up working directory
    """
    for file in glob.glob(os.path.join(work_path, '*')):
        move_status = False
        for save_type in save_types:
            if file.endswith(save_type):
                move_status = True
        if move_status:
            file_name = file.split('/')[-1]
            move_fpath = move_path+file_name
            shutil.move(file, move_fpath)
        else:
            os.remove(file)
    return


def execute(base_path, sites_path, work_path, library_path, recent_path, ddc_fpath, ddclist_fpath, run_combos):
    # determine all model runs and prompt user to proceed
    print "The following site / equilibrium / scenario combinations will be simulated: "
    for run_combo in run_combos:
        print "   %s / %s / %s" % (run_combo[0], run_combo[1], run_combo[2])
    proceed = raw_input("Select 'q' to quit or 'enter' to proceed: ")
    if proceed == 'q':
        return
    print

    # clear out current_results directory
    for file in glob.glob(os.path.join(recent_path, '*.*')):
        os.remove(file)

    # execute simulations
    start = time.time()
    for run_combo in run_combos:

        # read run specifications
        site = run_combo[1]
        bin_file = run_combo[2].split('.')[0]+'.bin'
        sch_file = run_combo[3]
        detail = run_combo[12]
        bin_handle = site+'_'+bin_file.split('.')[0]
        handle = bin_handle+'_'+sch_file.split('.')[0]
        print "Simulating %s / %s / %s..." % (site, bin_file, sch_file)
        print

        # define all input file paths
        schedule_fpath = base_path+sch_file
        site_path = sites_path+site+'/'
        site_fpath = site_path+bin_handle+'.100'
        weather_fpath = site_path+'weather.wth'
        binary_fpath = site_path+site+'_'+bin_file
        if not os.path.exists(binary_fpath):
            print "*** ERROR: No spin-up exists at ", binary_fpath
            print
            return

        # move files, execute simulations, and clean up
        copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, library_path)
        os.rename(work_path+bin_handle+'.100', work_path+'site.100')
        if detail:
            print "Running diagnostic model..."
            # DDcentEVI_diag_022715  # DDcentEVI_diag_112414'
            DDcentEVI('/data/paustian/GEVO_jeff_kent/modelDev/DDcentEVI_diag_022715',
                      sch_file, work_path, handle, ['dc_sip.csv', 'nflux.out'], ddclist_fpath=ddclist_fpath)
            DDC_detail(work_path, recent_path, handle, 2040, 2045)
        else:
            print "Running standard model..."
            DDcentEVI(ddc_fpath, sch_file, work_path, handle, ['dc_sip.csv', 'nflux.out'], ddclist_fpath=ddclist_fpath)
        clean_out(work_path, recent_path, ['.bin', '.lis', 'dc_sip.csv', 'nflux.out'])
        print
        print

    # report out
    print "****************************************************************************"
    print
    time_sec = round((time.time() - start), 2)
    time_min = round(time_sec/60.0, 2)
    print "It took %.2f minutes total to execute the specified DayCent simulation ensemble." % time_min
    print
    return


def spinup(base_path, sites_path, work_path, library_path, ddc_fpath, ddclist_fpath, run_combos):
    """
    """
    # determine spin-ups to run, and associated input files
    spinup_combos = []
    for run_combo in run_combos:
        spinup_combos.append((run_combo[1], run_combo[2]))
    unique_spinup_combos = list(set(spinup_combos))
    unique_spinup_combos.sort()
    print "Please specify a spin-up scenario to run, or 'a' for all:"
    for i, combo in enumerate(unique_spinup_combos):
        print "   ", i, combo[0], "-", combo[1]
    index = raw_input("")
    runs = []
    if index == 'a':
        runs = unique_spinup_combos
    else:
        index = int(index)
        runs.append(unique_spinup_combos[index])

    # execute all specified spin-ups
    for run in runs:
        site = run[0]
        sch_file = run[1]
        print "Running spin-up %s for site %s..." % (sch_file, site)

        # determine all necessary paths
        site_path = sites_path+site+'/'
        site_fpath = site_path+'site.100'
        weather_fpath = site_path+'weather.wth'
        schedule_fpath = base_path+sch_file

        # execute spin-up and plot som to verify equilibrium
        start = time.time()
        copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, library_path)
        handle = site+'_'+sch_file.split('.')[0]
        DDcentEVI(ddc_fpath, sch_file, work_path, handle, [], ddclist_fpath=ddclist_fpath, site_file_out=True)

        results = read_full_out(work_path+handle+'.lis', 3, 1)
        years = []
        for year in range(len(results[1])):
            years.append(year+1)
        plt.plot(years, list(results[1]))
        plt.xlabel('simulation time (years)')
        plt.ylabel('total SOM (gC/m2)')
        plt.savefig(handle+'.png')
        plt.close()

        # clean up and report out
        clean_out(work_path, site_path, ['.bin', '.lis', '.png', '%s.100' % handle])
        time_sec = round((time.time() - start), 2)
        time_min = round(time_sec/60.0, 2)
        print "It took %.2f minutes total to run the spin-up." % time_min
        print
        print
    return


def analysis(recent_path, archive_path, run_combos, labels):
    # create results archive
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H.%M")
    print "Please enter a short descriptive title for this analysis run, using underscores"
    description = raw_input("in place of spaces:  ")
    results_path = archive_path+timestamp+"_"+description+"/"
    os.mkdir(results_path)

    # copy over any growth visualization results from the current results directory
    for file in glob.glob(os.path.join(recent_path, '*')):
        if file.endswith('.png'):
            shutil.copy(file, results_path)

    # set up results output files
    out_fpath = results_path+"Table.csv"
    header = ["Simulation",
              "Avg. harvestable NPP (MgC/ha/y)",
              "Avg. SOC change (MgC/ha/y)",
              "Avg. system C change (MgCO2eq/ha/y)",
              "Avg. biomass harvest (Mg biomass/ha/y)",
              "Avg. bioenergy GHGE avoidance (MgCO2eq/ha/y)",
              "Avg. total GHGE mitigation (MgCO2eq/ha/y)"]
    out_file_object = open(out_fpath, "wb")
    c = csv.writer(out_file_object)
    c.writerow(header)

    bar_out_fpath = results_path+"Bar_plot.csv"
    bar_header = ["Avg. BGC accumulation (MgC/ha/y)",
                  "Avg. AGC accumulation (MgC/ha/y)",
                  "Avg. harvest (MgC/ha/y)",
                  "Avg. respiration & losses (MgC/ha/y)",
                  "Avg. NPP (MgC/ha/y)",
                  "Ecosystem C mitigation (MgCO2eq/ha/y)",
                  "Ecosystem N2O mitigation (MgCO2eq/ha/y)",
                  "Net bioenergy GHG mitigation (MgCO2eq/ha/y)",
                  "Supply chain GHG footprint (MgCO2eq/ha/y)",
                  "CCS marginal sequestration (MgCO2eq/ha/y)",
                  "Simulation"]
    bar_out_file_object = open(bar_out_fpath, "wb")
    d = csv.writer(bar_out_file_object)
    d.writerow(bar_header)

    raw_out_fpath = results_path+"Raw_plot.csv"
    raw_out_file_object = open(raw_out_fpath, "wb")
    e = csv.writer(raw_out_file_object)

    # initialize Figures
    fig0 = plt.figure(0)
    ax1 = fig0.add_subplot(2, 2, 1)
    ax2 = fig0.add_subplot(2, 2, 2, sharex=ax1)
    ax3 = fig0.add_subplot(2, 2, 3, sharex=ax1)

    fig1 = plt.figure(1)
    # gs1 = gridspec.GridSpec(2, 1)
    # gs1.update(right=0.4)
    # gs2 = gridspec.GridSpec(4, 1)
    # gs2.update(left=0.55, hspace=0.45)
    # # ax4 = plt.subplot(gs1[0, 0])
    # # ax5 = plt.subplot(gs2[0, 0])
    # # ax6 = plt.subplot(gs2[1, 0], sharex=ax5, sharey=ax5)
    ax7 = fig1.add_subplot(1, 1, 1)

    fig3 = plt.figure(3)
    ax8 = fig3.add_subplot(2, 4, 1)
    ax9 = fig3.add_subplot(2, 4, 2, sharex=ax8, sharey=ax8)
    ax10 = fig3.add_subplot(2, 4, 3, sharex=ax8, sharey=ax8)

    # define structure for saving stacked bar chart data [[avg_BGs], [avg_AGs], [avg_harvs], [avg_resp],
    #                                                     [ecosystem GHGE mitigation], [bioenergy GHGEM], [ticks]]
    bar_data = [[], [], [], [], [], [], [], [], [], [], []]

    # define structure for storing running curve minima and maxima for filling
    fill = {
            "saddlebrown": [ [[],[]], [[],[]], [[],[]] ],
            "g": [ [[],[]], [[],[]], [[],[]] ],
            "navy": [ [[],[]], [[],[]], [[],[]] ],
            "deepskyblue": [ [[],[]], [[],[]], [[],[]] ],
            "royalblue": [ [[],[]], [[],[]], [[],[]] ]
           }
    # panel12_fill = {
    #                 "m": True,
    #                 "g": True,
    #                 "r": False,
    #                 "b": True
    #                 }

    # define structure for saving raw daily mitigation data
    raw_data = []

    # step through run_combos, read associated output files, process results, and add to figures as appropriate
    counter = 0
    NPP_list = []
    for run_combo in run_combos:
        system, site, eq_file, sch_file, tech, panel12, panel5, panel6, panel3, color, linestyle, linewidth, detail,\
            panel8, panel9, panel10, color2, linestyle2,zorder = run_combo
        handle = site+'_'+eq_file.split('.')[0]+'_'+sch_file.split('.')[0]
        print "Processing case %s..." % handle
        lis_fpath = recent_path+handle+'.lis'
        nflux_fpath = recent_path+handle+'_nflux.out'
        dc_sip_fpath = recent_path+handle+'_dc_sip.csv'
        summary = [handle]   # data structure for tabular output

        ### process ANNUAL .lis output ################################################################################
        results = read_full_out(lis_fpath, 4, 1)
        years = []
        for year in range(len(results[0])):
            years.append(year+1)
        # time = np.array(results[0])
        somsc = np.array(results[1])
        crmvst = np.array(results[2])
        # aglivc = np.array(results[3])
        # stdedc = np.array(results[4])
        # fbrchc = np.array(results[5])
        # rlwodc = np.array(results[6])
        # wood1c = np.array(results[7])
        # wood2c = np.array(results[8])
        totsysc = np.array(results[9])
        # agcacc = np.array(results[10])
        # fbracc = np.array(results[11])
        # rlwacc = np.array(results[12])
        tcrem = np.array(results[13])

        # compute cumulative harvestable (aboveground) C production from ANNUAL output
        # if system == 'grass':
        #     cumulative_harvestable_gC_m2 = np.cumsum(agcacc)
        #     cumulative_harvestable_MgC_ha = cumulative_harvestable_gC_m2 * 0.01
        # elif system == 'forest':
        #     cumulative_harvestable_gC_m2 = fbrchc
        #     cumulative_harvestable_gC_m2 += rlwodc
        #     cumulative_harvestable_gC_m2 += wood1c
        #     cumulative_harvestable_gC_m2 += wood2c
        #     cumulative_harvestable_MgC_ha = cumulative_harvestable_gC_m2 * 0.01
        # elif system == 'bm':
        #     cumulative_harvestable_gC_m2 = np.cumsum(crmvst)
        #     cumulative_harvestable_MgC_ha = cumulative_harvestable_gC_m2 * 0.01
        # summary.append(round(cumulative_harvestable_MgC_ha[-1]/float(len(cumulative_harvestable_MgC_ha)), 3))

        # compute SOC and total ecosystem C from ANNUAL output
        tot_SOC_MgC_ha = somsc * 0.01
        tot_sysC_MgC_ha = totsysc * 0.01

        # compute harvest from ANNUAL output
        wood_harvest = np.array(tcrem[0])
        wood_harvest = np.append(wood_harvest, np.diff(tcrem))
        harvest_gc_m2 = (crmvst + wood_harvest)
        nonzero_harvest_gc_m2 = []
        for harvest in harvest_gc_m2:
            if harvest > 0:
                nonzero_harvest_gc_m2.append(harvest)
        avg_harvest_gC_m2 = np.sum(nonzero_harvest_gc_m2) / float(len(nonzero_harvest_gc_m2))

        ### process DAILY dc_sip.csv output ###########################################################################
        more_results = read_full_out(dc_sip_fpath, 1, 0, delim='c')
        days = []
        zeros = []
        for day in range(len(more_results[0])):
            days.append((day+1)/365.0)
            zeros.append(0)
        mcprd1 = np.array(more_results[22])   # Daily NPP for shoots for grass/crop system (gC/m2)
        mfprd1 = np.array(more_results[25])   # Daily NPP for live leaves for tree system (gC/m2)
        mfprd3 = np.array(more_results[28])   # Daily NPP for live fine branches for tree system (gC/m2)
        mfprd4 = np.array(more_results[29])   # Daily NPP for live large wood for tree system (gC/m2)
        NPP = np.array(more_results[31])
        aglivc = np.array(more_results[33])
        bglivcj = np.array(more_results[34])
        bglivcm = np.array(more_results[35])
        rleavc = np.array(more_results[36])
        frootcj = np.array(more_results[37])
        frootcm = np.array(more_results[38])
        fbrchc = np.array(more_results[39])
        rlwodc = np.array(more_results[40])
        crootc = np.array(more_results[41])
        tlai = np.array(more_results[42])
        stdedc = np.array(more_results[43])
        wood1c = np.array(more_results[44])   # dead fine branches
        wood2c = np.array(more_results[45])   # dead large wood
        wood3c = np.array(more_results[46])   # dead coarse roots
        strucc1 = np.array(more_results[47])  # 1=surface, 2=soil
        metabc1 = np.array(more_results[48])
        strucc2 = np.array(more_results[49])
        metabc2 = np.array(more_results[50])
        som1c1 = np.array(more_results[51])
        som1c2 = np.array(more_results[52])
        som2c1 = np.array(more_results[53])
        som2c2 = np.array(more_results[54])
        som3c = np.array(more_results[55])

        # compute aboveground and belowground totals from DAILY output
        AG = aglivc
        AG += rleavc
        AG += fbrchc
        AG += rlwodc
        AG += stdedc
        AG += wood1c
        AG += wood2c
        AG += strucc1
        AG += metabc1

        BG = bglivcj
        BG += bglivcm
        BG += frootcj
        BG += frootcm
        BG += crootc
        BG += wood3c
        BG += strucc2
        BG += metabc2
        BG += som1c1
        BG += som1c2
        BG += som2c1
        BG += som2c2
        BG += som3c

        # compute ANPP from DAILY output
        ANPP = mcprd1
        ANPP += mfprd1
        ANPP += mfprd3
        ANPP += mfprd4

        # compute cumulative sum arrays, net change arrays, and perform unit conversions from DAILY output
        cum_NPP = np.cumsum(NPP) * 0.01
        cum_ANPP = np.cumsum(ANPP) * 0.01
        initial_AG = AG[0]
        net_AG = np.array([x-initial_AG for x in AG]) * 0.01
        initial_BG = BG[0]
        net_BG = np.array([x-initial_BG for x in BG]) * 0.01
        cum_loss = cum_NPP - net_AG
        cum_loss -= net_BG
        net_eco_MgC_ha = net_AG + net_BG
        net_eco_MgC_ha.tolist()
        eco_C_MgCO2_ha = [x*3.67 for x in net_eco_MgC_ha]
        avg_ecosystem_C_mitigation = eco_C_MgCO2_ha[-14975] / float(days[-14975])   # i.e., including only first 30 y of simulation, skipping last 41
        bar_data[5].append(avg_ecosystem_C_mitigation)

        summary.append(round(cum_ANPP[-1] / float(days[-1]), 3))
        summary.append(round(np.mean(np.diff(tot_SOC_MgC_ha)), 4))
        summary.append(round(np.mean(np.diff(tot_sysC_MgC_ha)) * 3.67, 3))
        summary.append(round((avg_harvest_gC_m2*0.01)/0.45, 3))

        ghge_MgCO2_ha = np.array(eco_C_MgCO2_ha)

        ### process DAILY nflux.out output ####################################################################
        n2o_results = read_full_out(nflux_fpath, 1, 0)
        nitri = np.array(n2o_results[2])   # Daily total nitrification N2O emissions (gN/ha) not /m2!!!
        denitri = np.array(n2o_results[3])   # Daily total denitrification N2O emissions (gN/ha) not /m2!!!
        N2Oflux = nitri
        N2Oflux += denitri
        N2O_MgCO2_ha = N2Oflux * 0.000001 * n_to_n2o * n2o_gwp100 * -1.0
        cum_N2O_MgCO2_ha = np.cumsum(N2O_MgCO2_ha)
        avg_ecosystem_N_mitigation = cum_N2O_MgCO2_ha[-14975] / float(days[-14975])
        bar_data[6].append(avg_ecosystem_N_mitigation)

        # print len(ghge_MgCO2eq_ha)
        # print len(cum_N2O_MgCO2eq_ha)
        ghge_MgCO2_ha += cum_N2O_MgCO2_ha

        ### compute and plot mitigation from ecosystem C changes and bioenergy production #############################
        # update harvest array from annual to daily basis
        daily_harvest_g_m2 = []
        daily_counter = 0
        annual_counter = 0
        harvest_doy = 230
        for k in range(len(days)):
            if k % 365 == 0:
                daily_counter = 0
            if daily_counter == harvest_doy:
                daily_harvest_g_m2.append(harvest_gc_m2[annual_counter])
                annual_counter += 1
            else:
                daily_harvest_g_m2.append(0)
            daily_counter += 1
        daily_harvest_g_m2 = np.array(daily_harvest_g_m2)

        # compute bioenergy GHG mitigation
        if tech:
            # if system == "forest":
            #     print "   harvest", harvest_gc_m2
            logistic_efficiency = 1.0   # previously 0.9
            daily_harvest_g_m2 = daily_harvest_g_m2 * logistic_efficiency
            net_ghg_avoid_gCO2_m2, supply_chain_gCO2_m2, CCS_marg_gCO2_m2 = bioenergy_ghg_avoid_gCO2eq_m2(daily_harvest_g_m2, tech)

            nonzero_net_ghg_avoid_gCO2_m2 = []
            for ghge in net_ghg_avoid_gCO2_m2:
                if ghge > 0:
                    nonzero_net_ghg_avoid_gCO2_m2.append(ghge)
            summary.append(round((np.sum(nonzero_net_ghg_avoid_gCO2_m2)/float(len(nonzero_net_ghg_avoid_gCO2_m2)))*0.01, 3))

            cum_net_ghg_avoid_gCO2_m2 = np.cumsum(net_ghg_avoid_gCO2_m2)
            cum_supply_chain_gCO2_m2 = np.cumsum(supply_chain_gCO2_m2)
            cum_CCS_marg_gCO2_m2 = np.cumsum(CCS_marg_gCO2_m2)

            cum_ghg_avoid_MgCO2_ha = cum_net_ghg_avoid_gCO2_m2 * 0.01
            cum_supply_chain_MgCO2_ha = cum_supply_chain_gCO2_m2 * 0.01
            cum_CCS_marg_MgCO2_ha = cum_CCS_marg_gCO2_m2 * 0.01

            ghge_MgCO2_ha += cum_ghg_avoid_MgCO2_ha
            avg_bioenergy_mitigation = cum_ghg_avoid_MgCO2_ha[-14975] / float(days[-14975])
            avg_supply_chain = cum_supply_chain_MgCO2_ha[-14975] / float(days[-14975])
            avg_CCS_marg = cum_CCS_marg_MgCO2_ha[-14975] / float(days[-14975])
            bar_data[7].append(avg_bioenergy_mitigation)
            bar_data[8].append(avg_supply_chain)
            bar_data[9].append(avg_CCS_marg)
        else:
            summary.append(0)
            bar_data[7].append(0)
            bar_data[8].append(0)
            bar_data[9].append(0)

        summary.append(round(ghge_MgCO2_ha[-1] / float(days[-1]), 3))

        # record raw daily GHG mitigation output to file
        if not counter:
            raw_data.append(['Years post-conversion']+days)
        tabular_ghge_MgCO2_ha = ghge_MgCO2_ha.tolist()
        raw_data.append([handle]+tabular_ghge_MgCO2_ha)

        # record results to file, and add to plot where appropriate
        # truncate bioenergy results for improved read-ability
        if tech:
            ghge_MgCO2_ha_plot = ghge_MgCO2_ha[:-1200]
            days_plot = days[:-1200]
        else:
            ghge_MgCO2_ha_plot = ghge_MgCO2_ha
            days_plot = days
        c.writerow(summary)
        plt.figure(0)
        if panel12:
            ax1.plot(days, cum_ANPP, color=color, linestyle=linestyle)
            # if list(fill[color][0][0]):
            #     fill[color][0][0] = np.minimum(fill[color][0][0], cum_ANPP)
            # else:
            #     fill[color][0][0] = cum_ANPP
            # if list(fill[color][0][1]):
            #     fill[color][0][1] = np.maximum(fill[color][0][1], cum_ANPP)
            # else:
            #     fill[color][0][1] = cum_ANPP
            ax2.plot(years, tot_sysC_MgC_ha, color=color, linestyle=linestyle)
            # if list(fill[color][1][0]):
            #     fill[color][1][0] = np.minimum(fill[color][1][0], tot_sysC_MgC_ha)
            # else:
            #     fill[color][1][0] = tot_sysC_MgC_ha
            # if list(fill[color][1][1]):
            #     fill[color][1][1] = np.maximum(fill[color][1][1], tot_sysC_MgC_ha)
            # else:
            #     fill[color][1][1] = tot_sysC_MgC_ha
        if panel3:
            ax3.plot(days, ghge_MgCO2_ha, color=color, linestyle=linestyle)
            ax7.plot(days, ghge_MgCO2_ha, color=color, linestyle=linestyle)
            # if list(fill[color][2][0]):
            #     fill[color][2][0] = np.minimum(fill[color][2][0], ghge_MgCO2eq_ha)
            # else:
            #     fill[color][2][0] = ghge_MgCO2eq_ha
            # if list(fill[color][2][1]):
            #     fill[color][2][1] = np.maximum(fill[color][2][1], ghge_MgCO2eq_ha)
            # else:
            #     fill[color][2][1] = ghge_MgCO2eq_ha

        if panel8:
            ax8.plot(days_plot, ghge_MgCO2_ha_plot, color=color2, linestyle=linestyle2, linewidth=1, zorder=zorder)
            if list(fill[color2][0][0]):
                fill[color2][0][0] = np.minimum(fill[color][0][0], ghge_MgCO2_ha_plot)
            else:
                fill[color2][0][0] = ghge_MgCO2_ha_plot
            if list(fill[color2][0][1]):
                fill[color2][0][1] = np.maximum(fill[color][0][1], ghge_MgCO2_ha_plot)
            else:
                fill[color2][0][1] = ghge_MgCO2_ha_plot
        if panel9:
            ax9.plot(days_plot, ghge_MgCO2_ha_plot, color=color2, linestyle=linestyle2, linewidth=1, zorder=zorder)
            if list(fill[color2][1][0]):
                fill[color2][1][0] = np.minimum(fill[color][1][0], ghge_MgCO2_ha_plot)
            else:
                fill[color2][1][0] = ghge_MgCO2_ha_plot
            if list(fill[color2][1][1]):
                fill[color2][1][1] = np.maximum(fill[color][1][1], ghge_MgCO2_ha_plot)
            else:
                fill[color2][1][1] = ghge_MgCO2_ha_plot
        if panel10:
            ax10.plot(days_plot, ghge_MgCO2_ha_plot, color=color2, linestyle=linestyle2, linewidth=1, zorder=zorder)
            if list(fill[color2][2][0]):
                fill[color2][2][0] = np.minimum(fill[color][2][0], ghge_MgCO2_ha_plot)
            else:
                fill[color2][2][0] = ghge_MgCO2_ha_plot
            if list(fill[color2][2][1]):
                fill[color2][2][1] = np.maximum(fill[color][2][1], ghge_MgCO2_ha_plot)
            else:
                fill[color2][2][1] = ghge_MgCO2_ha_plot

        # add to stack bar chart data structure
        avg_NPP = cum_NPP[-14975] / float(days[-14975])
        avg_BG = net_BG[-14975] / float(days[-14975])
        avg_AG = net_AG[-14975] / float(days[-14975])
        avg_harv = np.mean(harvest_gc_m2[:14975]) * 0.01
        avg_resp = (cum_loss[-14975] / float(days[-14975])) - avg_harv
        bar_data[0].append(avg_BG)
        bar_data[1].append(avg_AG)
        bar_data[2].append(avg_harv)
        bar_data[3].append(avg_resp)
        bar_data[4].append(avg_NPP)
        NPP_list.append(avg_NPP)
        bar_data[10].append(handle)

        counter += 1

        # plot dynamics details where appropriate
        start = 0
        end = 1200
        # if panel5:
        #     stacked_line_plot(ax5,
        #                       days[start:end],
        #                       [zeros[start:end], net_BG[start:end], net_AG[start:end], cum_loss[start:end]],
        #                       [('g', '', 'BGC change'), ('g', '///', 'AGC change'), ('w', '', 'respiration/losses')],
        #                       legend=True)
        # if panel6:
        #     stacked_line_plot(ax6,
        #                       days[start:end],
        #                       [zeros[start:end], net_BG[start:end], net_AG[start:end], cum_loss[start:end]],
        #                       [('m', '', 'BGC change'), ('m', '///', 'AGC change'), ('w', '', 'respiration/losses')])

    ### finalize figure formatting ###################################################################################
    # create upper left stacked bar chart
    bar_colors_labels = [('k', 'belowground C'), ('c', 'aboveground C'), ('r', 'harvest'), ('w', 'respiration/losses')]
    increment = range(len(bar_data[0]))
    width = 0.75
    bottoms = []
    for a in increment:
        bottoms.append(0)
    for i, color_label in enumerate(bar_colors_labels):
        color, label = color_label
        # ax4.bar(increment, bar_data[i], width=width, bottom=bottoms, color=color, label=label)
        for j in range(len(bottoms)):
            bottoms[j] += bar_data[i][j]
    for k in increment:
        increment[k] += width/2.0
    # ax4.set_xticks(increment, bar_data[4])
    # ax4.set_xticklabels(bar_data[4], rotation=90, )
    # ax4.legend(loc=2, prop={'size': 7.5})

    # write bar chart output to file for use in Excel
    bar_data = zip(*bar_data)
    for line in bar_data:
        d.writerow(line)

    # write bar chart output to file for use in Excel
    raw_data = zip(*raw_data)
    for line in raw_data:
        e.writerow(line)

    # add fill to lower left time-series plot
    # print
    for color in fill:
        # print "Filling color: ", color
        if list(fill[color][0][0]):
            # print "   Filling ax8"
            # print "      Min = ", fill[color][0][0][-1]
            # print "      Max = ", fill[color][0][1][-1]
            days = []
            for day in range(len(fill[color][0][0])):
                days.append((day+1)/365.0)
            ax8.fill_between(days, fill[color][0][0], fill[color][0][1], facecolor=color, alpha=0.5, linewidth=0.0, zorder=0)
        if list(fill[color][1][0]):
            # print "   Filling ax9"
            # print "      Min = ", fill[color][1][0][-1]
            # print "      Max = ", fill[color][1][1][-1]
            days = []
            for day in range(len(fill[color][1][0])):
                days.append((day+1)/365.0)
            ax9.fill_between(days, fill[color][1][0], fill[color][1][1], facecolor=color, alpha=0.5, linewidth=0.0, zorder=0)
        if list(fill[color][2][0]):
            # print "   Filling ax10"
            # print "      Min = ", fill[color][2][0][-1]
            # print "      Max = ", fill[color][2][1][-1]
            days = []
            for day in range(len(fill[color][2][0])):
                days.append((day+1)/365.0)
            ax10.fill_between(days, fill[color][2][0], fill[color][2][1], facecolor=color, alpha=0.5, linewidth=0.0, zorder=0)
            # ax7.fill_between(days, fill[color][2][0], fill[color][2][1], facecolor=color, alpha=0.15, linewidth=0.0)
    # print

    # update plot formatting and add axis labels
    matplotlib.rcParams.update({'font.size': 12})
    ax1.set_title("Cumulative ANPP", size=11, fontweight='bold')
    ax1.set_ylabel("MgC / ha")
    ax2.set_title("Current total ecosystem carbon", size=11, fontweight='bold')
    ax2.set_ylabel("MgC / ha")
    ax3.set_title("Cumulative GHGE benefit", size=11, fontweight='bold')
    ax3.set_ylabel("MgCO2eq / ha")
    ax3.set_xlabel("Time after harvest/LUC (years)", size=10)
    # ax4.set_title("Avgeraged annual carbon fluxes",  size=11, fontweight='bold')
    # ax4.set_ylabel("MgC/ha/y", size=10)
    # ax5.set_title("Carbon flux detail",  size=11, fontweight='bold')
    # ax5.set_ylabel("MgC/ha                                  .", size=10)
    # ax5.set_xlabel("Time after harvest/LUC (years)", size=10)
    # ax5.set_ylim(0, 20)
    # ax7.set_title("GHGE mitigation detail", size=11, fontweight='bold')
    ax7.set_ylabel('GHGE (Mg $\mathregular{CO_{2}}$eq $\mathregular{(ha)^{-1}}$)', size=12)
    ax7.set_xlabel("Time after harvest or land use change (years)", size=12)
    # ax5.text(0.42, 0.8, "Unmanaged grassland", transform=ax5.transAxes, fontsize=9, fontweight='bold')
    # ax6.text(0.04, 0.8, "Unmanaged forest", transform=ax6.transAxes, fontsize=9, fontweight='bold')
    ax1.axhline(color='k')
    ax2.axhline(color='k')
    ax3.axhline(color='k')
    # ax4.axhline(color='k')
    # ax5.axhline(color='k')
    # ax6.axhline(color='k')
    ax7.axhline(color='k')
    ax7.set_ylim(-1000, 1500)

    # add timeseries subplot titles and horizontal lines
    ax8.axhline(color='k', zorder=5)
    ax8.set_xlabel("Years since conversion")
    ax8.set_ylabel("GHGE mitigation (Mg $\mathregular{CO_{2}}$eq $\mathregular{(ha)^{-1}}$)")
    ax8.set_title("Crop/pasture,\notherwise forest")
    # ax8.set_title("Crop or pasture\nto biofuel")
    ax9.axhline(color='k', zorder=5)
    ax9.set_title("Initial Land Use:\nCrop/pasture,\notherwise grassland")
    # ax9.set_title("Crop or pasture\nto regrowth")
    ax10.axhline(color='k', zorder=5)
    ax10.set_title("70 y.o. forest")
    # ax10.set_title("70 y.o. forest\nto biofuel or\nregrowth")

    # timeseries subplot formatting - declutter redundant axes, but add grids
    ax8.spines['right'].set_color('none')
    ax8.spines['top'].set_color('none')
    ax8.get_xaxis().tick_bottom()
    ax8.get_yaxis().tick_left()
    ax8.xaxis.set_ticks([0, 20, 40, 60, 80])
    ax8.grid(True)

    # ax9.set_frame_on(False)
    ax9.spines['left'].set_color('none')
    ax9.spines['right'].set_color('none')
    ax9.spines['top'].set_color('none')
    ax9.get_xaxis().tick_bottom()
    ax9.get_yaxis().tick_left()
    ax9.grid(True)
    # plt.setp(ax9.get_xticklabels(), visible=False)
    plt.setp(ax9.get_yticklabels(), visible=False)

    # ax10.set_frame_on(False)
    ax10.spines['left'].set_color('none')
    ax10.spines['right'].set_color('none')
    ax10.spines['top'].set_color('none')
    ax10.get_xaxis().tick_bottom()
    ax10.get_yaxis().tick_left()
    ax10.grid(True)
    # plt.setp(ax10.get_xticklabels(), visible=False)
    plt.setp(ax10.get_yticklabels(), visible=False)

    # add legend
    dummy_point = (-1, -1)
    for entry in labels:
        label, color, linestyle, linewidth = entry
        ax3.plot(dummy_point[0], dummy_point[1], label=label, color=color, linestyle=linestyle)
    for entry in labels:
        label, color, linestyle, linewidth = entry
        ax7.plot(dummy_point[0], dummy_point[1], label=label, color=color, linestyle=linestyle, linewidth=linewidth)
        ax9.plot(dummy_point[0], dummy_point[1], label=label, color=color, linestyle=linestyle, linewidth=10)
    ax3.legend(bbox_to_anchor=(1.05, 0.9), loc=2, prop={'size': 8.5})
    # ax7.legend(loc=4, prop={'size': 6.5})
    ax7.legend(loc=2, prop={'size': 12})
    ax9.legend(bbox_to_anchor=(0.5, -0.9), loc=8, prop={'size': 9}, title="Subsequent Land Use:")

    # save figure, table, and report out
    out_file_object.close()
    shutil.copy(results_path+'Table.csv', archive_path+'Table.csv')
    bar_out_file_object.close()
    shutil.copy(results_path+'Bar_plot.csv', archive_path+'Bar_plot.csv')
    raw_out_file_object.close()
    shutil.copy(results_path+'Raw_plot.csv', archive_path+'Raw_plot.csv')
    plt.figure(0)
    plt.savefig(results_path+'Fig3.png')
    shutil.copy(results_path+'Fig3.png', archive_path+'Fig3.png')
    plt.close()
    plt.figure(1)
    plt.savefig(results_path+'Fig3_alternate.png')
    shutil.copy(results_path+'Fig3_alternate.png', archive_path+'Fig3_alternate.png')
    plt.close()
    plt.figure(3)
    matplotlib.rcParams.update({'font.size': 9.5})
    plt.savefig(results_path+'Fig2B-D.png')
    shutil.copy(results_path+'Fig2B-D.png', archive_path+'Fig2B-D.png')
    plt.close()
    results_path = results_path.split('/')[-2]
    print
    print "Resulting figure and tabular output archived in the directory "+results_path
    print
    # print "Scenario average NPP values (Mg C ha-1 y-1)"
    # print NPP_list
    # print "   Last/first ratio: %.3f" % (NPP_list[-1]/float(NPP_list[0]))
    # print

    return


def calibration(base_path, cal_sets, calibration_path, work_path, library_path, ddc_fpath, ddclist_fpath):
    # determine all calibration runs and prompt user to proceed
    print "Calibration will be performed for the following regions: "
    cal_paths = filter(os.path.isdir, [os.path.join(calibration_path, f) for f in os.listdir(calibration_path)])
    cal_regions = []
    for cal_path in cal_paths:
        cal_run = cal_path.split('/')[-1]
        print "   %s" % cal_run
        cal_regions.append(cal_run)
    proceed = raw_input("Select 'q' to quit, 's' to re-run the calibration spin-ups, or 'enter' to proceed: ")
    if proceed == 'q':
        return
    if proceed == 's':
        # determine spin-ups to run, and associated input files
        print
        print "Please specify the number for a specific region to spin up, or 'a' for all:"
        for i, combo in enumerate(cal_paths):
            print "   ", i, combo
        index = raw_input("")
        regions = []
        if index == 'a':
            regions = cal_paths
        else:
            print index
            index = int(index)
            print type(index)
            print cal_paths[index]
            regions.append(cal_paths[index])

        # execute spin-ups for the specified regions
        for i, cal_path in enumerate(regions):
            sch_file = cal_path.split('/')[-1]+'_pine_spinup.sch'
            print "Executing spin-ups for %s" % cal_path
            print

            # determine all necessary paths
            inter_path = cal_path+'/'
            site_paths = filter(os.path.isdir, [os.path.join(inter_path, f) for f in os.listdir(inter_path)])
            print "Spinning up the following sites:"
            for site_path in site_paths:
                print site_path
            print

            for site_path in site_paths:

                site_path += '/'
                site_fpath = site_path+'site.100'
                weather_fpath = site_path+'weather.wth'
                schedule_fpath = base_path+sch_file

                # execute spin-up and plot som to verify equilibrium
                start = time.time()
                copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, library_path)
                handle = 'spinup'
                DDcentEVI(ddc_fpath, sch_file, work_path, handle, [], ddclist_fpath=ddclist_fpath, site_file_out=True)

                results = read_full_out(work_path+handle+'.lis', 3, 1)
                years = []
                for year in range(len(results[1])):
                    years.append(year+1)
                plt.plot(years, list(results[1]))
                plt.xlabel('simulation time (years)')
                plt.ylabel('total SOM (gC/m2)')
                plt.savefig(handle+'.png')
                plt.close()

                # clean up and report out
                clean_out(work_path, site_path, ['.bin', '.lis', '.png', '%s.100' % handle])
                time_sec = round((time.time() - start), 2)
                time_min = round(time_sec/60.0, 2)
                print "It took %.2f minutes total to run the spin-up." % time_min
                print
                print
    print

    # clear out current_results directory
    for file in glob.glob(os.path.join(recent_path, '*.*')):
        os.remove(file)

    # prompt for archive label
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H.%M")
    print "Please enter a short descriptive title for this calibration run, using underscores"
    description = raw_input("in place of spaces:  ")

    # execute calibration simulations
    start = time.time()
    for cal_set in cal_sets:
        region, reference_values = cal_set
        sch_file = region+'_greenbook.sch'
        region_path = calibration_path+region+'/'
        schedule_fpath = base_path + '/' + sch_file
        cal_site_paths = filter(os.path.isdir, [os.path.join(region_path, f) for f in os.listdir(region_path)])

        for cal_site_path in cal_site_paths:
            print "Running calibration simulation %s..." % cal_site_path
            print

            # define all input file paths
            index = cal_site_path[-1:]
            cal_site_path += '/'
            site_fpath = cal_site_path+'spinup.100'
            weather_fpath = cal_site_path+'weather.wth'
            binary_fpath = cal_site_path+'spinup.bin'
            if not os.path.exists(binary_fpath):
                print "*** ERROR: No spin-up exists at ", binary_fpath
                print
                return

            # move files, execute simulations, and clean up
            copy_in(work_path, schedule_fpath, site_fpath, weather_fpath, library_path)
            os.rename(work_path+'spinup.100', work_path+'site.100')
            handle = region+'_'+index+'_'+sch_file[:-4]
            DDcentEVI(ddc_fpath, sch_file, work_path, handle, ['dc_sip.csv'], ddclist_fpath=ddclist_fpath)
            clean_out(work_path, recent_path, ['.bin', '.lis', 'dc_sip.csv'])
            print
            print

    # report out
    print "****************************************************************************"
    print
    time_sec = round((time.time() - start), 2)
    time_min = round(time_sec/60.0, 2)
    print "It took %.2f minutes total to execute the specified DayCent calibration ensemble." % time_min
    print

    # create results archive
    results_path = archive_path+timestamp+"_Calibration-"+description+'/'
    os.mkdir(results_path)
    for file in glob.glob(os.path.join(recent_path, '*')):
        file_name = file.split('/')[-1]
        move_fpath = results_path+file_name
        shutil.move(file, move_fpath)

    # step through calibration sets, read associated output files, process results, and add to figures as appropriate
    for cal_set in cal_sets:
        region, reference_values = cal_set
        deadc_Mg_ha_set = []
        totalc_Mg_ha_set = []
        region_path = calibration_path+region+'/'
        cal_site_paths = filter(os.path.isdir, [os.path.join(region_path, f) for f in os.listdir(region_path)])

        for cal_site_path in cal_site_paths:
            index = cal_site_path[-1:]
            cal_site_path += '/'
            handle = region+'_'+index+'_'+region+'_greenbook'
            print "Processing case %s..." % cal_site_path
            lis_fpath = results_path+handle+'.lis'
            results = read_full_out(lis_fpath, 4, 1)

            # annual output data processing
            years = []
            for year in range(len(results[0])):
                years.append(year+1)
            # time = np.array(results[0])
            # somsc = np.array(results[1])
            # crmvst = np.array(results[2])
            # aglivc = np.array(results[3])
            # stdedc = np.array(results[4])
            # fbrchc = np.array(results[5])
            rlwodc = np.array(results[6])
            wood1c = np.array(results[7])
            wood2c = np.array(results[8])
            # totsysc = np.array(results[9])
            # agcacc = np.array(results[10])
            fbracc = np.array(results[11])
            # rlwacc = np.array(results[12])
            # tcrem = np.array(results[13])
            rleavc = np.array(results[14])
            crootc = np.array(results[15])
            # wood3c = np.array(results[16])

            # compute live and dead biomass totals, units conversion
            deadc = wood1c   # g / m2
            deadc += wood2c

            totalc = deadc.copy()   # g / m2
            totalc += rlwodc
            totalc += fbracc
            totalc += rleavc
            totalc += crootc

            deadc_Mg_ha = deadc * 0.01
            totalc_Mg_ha = totalc * 0.01
            deadc_Mg_ha_set.append(deadc_Mg_ha)
            totalc_Mg_ha_set.append(totalc_Mg_ha)

        # report the spread in results
        print
        print "End-of-simulation variability in total ecosystem carbon (Mg C ha-1):"
        for totalc_Mg_ha in totalc_Mg_ha_set:
            print totalc_Mg_ha[-1]
        print

        # compute average simulated across sites and plot versus reference
        length = len(totalc_Mg_ha_set)
        deadc_Mg_ha_avg = 0
        totalc_Mg_ha_avg = 0
        for i in range(length):
            if i==0:
                deadc_Mg_ha_avg = deadc_Mg_ha_set[0]
                totalc_Mg_ha_avg = totalc_Mg_ha_set[0]
            else:
                deadc_Mg_ha_avg += deadc_Mg_ha_set[i]
                totalc_Mg_ha_avg += totalc_Mg_ha_set[i]
        deadc_Mg_ha_avg /= float(length)
        totalc_Mg_ha_avg /= float(length)

        # plot observed & modeled ecosystem carbon trends
        plt.plot(years, deadc_Mg_ha_avg, marker='None', linestyle='-', color='r', label='DayCent simulated dead C')
        plt.plot(years, totalc_Mg_ha_avg, marker='None', linestyle='-', color='g', label='DayCent simulated total (live+dead) C')

        ref_years = np.array(reference_values[0])
        ref_livec = np.array(reference_values[1])
        ref_deadc = np.array(reference_values[2])
        ref_totalc = ref_deadc.copy()
        ref_totalc += ref_livec
        plt.scatter(ref_years, ref_deadc, marker='o', color='r', label='Greenbook estimated dead C')
        plt.scatter(ref_years, ref_totalc, marker='o', color='g', label='Greenbook estimated total (live+dead) C')

        plt.xlabel('Stand age (y)')
        plt.ylabel('Stand biomass (MgC/ha)')
        plt.title('%s regional calibration' % region)
        plt.axhline(0, color='k')
        plt.axvline(70, color='grey', linewidth=10, alpha=0.3)
        plt.legend(loc=2, frameon=False, prop={'size': 11})
        plt.savefig(results_path+'%s_calibration.png' % region)
        plt.close()

        # report final productivity for calibration purposes
        print "Observed and modeled end-of-simulation productivity (Mg C ha-1)"
        print "   Observed: %.1f  Modeled: %.1f" % (ref_totalc[-1], totalc_Mg_ha_avg[-1])
        print "   Modeled/observed ratio: %.3f" % (totalc_Mg_ha_avg[-1]/float(ref_totalc[-1]))
        print


def soil_comparison(mukey_fpath, CFARM_db_fpath, base_path, archive_path, work_path, library_path, ddc_fpath, ddclist_fpath):

    # specify comparisons to be made: scenario title, eq.sch, and test.sch
    comparisons = [
        ('Unmanaged_pine_old_cal', 'pine_spinup.sch', 'greenbook.sch'),
        # ('Unmanaged_pine_new_cal', 'pine_spinup_DUKE2.sch', 'greenbook_DUKE2.sch'),
        ('Unmanaged_grass', 'grass_lowSOCcrop_eq.sch', 'grass_NY.sch',),
        ('Managed_switchrass', 'grass_lowSOCcrop_eq.sch', 'switchgrass_NY.sch',)
    ]

    # comparisons = [
    #     ('Unmanaged_pine_old_cal', 'greenbook.sch', 'greenbook2.sch'),
    #     ('Unmanaged_pine_new_cal', 'greenbook.sch', 'greenbook2.sch')
    # ]

    comparison_results = [
        [],
        [],
        [],
        []
    ]

    header = ['scenario', 'equilibrium_sch', 'test_sch']

    # specify paths
    wth_base_path = '/data/paustian/AFRI/NARR_gridxy_wth/'
    ddcarch_fpath = '/data/paustian/AFRI/simulations/run_ddc/model/DDCsitarchive'
    site_arch_fpath = '/data/paustian/AFRI/simulations/site100/runfile_cometfarm_spinup_current.dsa'
    fips = 36065
    wth_fpath = wth_base_path + "NARR_135_251.wth"   # previously "NARR_133_286.wth"

    # read and display all soils to be tested, prompt user to proceed, and create results archive
    print "The following comparisons will be made: "
    for comparison in comparisons:
        print "   %s (%s, %s)" % (comparison[0], comparison[1], comparison[2])
    print

    print "Comparison will be performed for the following soils: "
    mukey_fobj = open(mukey_fpath, 'rU')
    metadata = csv.reader(mukey_fobj)
    next(metadata)
    for row in metadata:
        print "   ", int(row[0])
    print
    description = raw_input("Please enter a descriptive title for the analysis: ")
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H.%M")
    results_path = archive_path+timestamp+"_Soil_productivity-"+description+'/'
    os.mkdir(results_path)

    # establish connection with CFARM archive database
    con = sqlite3.connect(CFARM_db_fpath)
    cur = con.cursor()

    # cycle through each soil
    mukey_fobj = open(mukey_fpath, 'rU')
    metadata = csv.reader(mukey_fobj)
    next(metadata)
    for row in metadata:
        mukey = int(row[0])
        header.append(mukey)
        print "*******************************************************************************************************"
        print "Analyzing soil: ", mukey

        # select a single representative CFARM archive runno, and write associated site and weather files to workspace
        cur.execute("""SELECT runno FROM CFARM_archive WHERE mukey=%i AND fips=%i AND irrig='N' AND grazetill='graze'
                       LIMIT 1""" % (mukey, fips))
        runno = cur.fetchall()[0][0]
        site_fpath = work_path+"site.100"
        site_command = "%s -r %s %s %s" % (ddcarch_fpath, runno, site_arch_fpath, site_fpath)

        # cycle through each comparison scenario
        for i, comparison in enumerate(comparisons):
            scenario, eq_schfile, test_schfile = comparison
            print "Executing spin-up..."
            print

            # copy over schedule file and library
            subprocess.call("%s" % site_command, shell=True)
            shutil.copy(wth_fpath, work_path+'weather.wth')
            schedule_fpath = base_path + eq_schfile
            shutil.copy(schedule_fpath, work_path)
            for file in glob.glob(os.path.join(library_path, '*')):
                shutil.copy(file, work_path)

            # execute spinup, plot SOC trend
            handle = scenario + '_' + str(mukey) + 'eq'
            DDcentEVI(ddc_fpath, schedule_fpath, work_path, handle, [], ddclist_fpath=ddclist_fpath, site_file_out=True)
            results = read_full_out(work_path+handle+'.lis', 3, 1)
            years = []
            for year in range(len(results[1])):
                years.append(year+1)
            plt.plot(years, list(results[1]))
            plt.xlabel('simulation time (years)')
            plt.ylabel('total SOM (gC/m2)')
            plt.savefig(results_path+handle+'.png')
            plt.close()

            # execute future simulation and record results
            print "Executing forward simulation..."
            print
            schedule_fpath = base_path + test_schfile
            shutil.copy(schedule_fpath, work_path)
            handle = scenario + '_' + str(mukey)
            DDcentEVI(ddc_fpath, schedule_fpath, work_path, handle, [], ddclist_fpath=ddclist_fpath)

            # analyze results & clean up
            print "Processing daily NPP results..."
            dc_sip_fpath = work_path + 'dc_sip.csv'
            more_results = read_full_out(dc_sip_fpath, 1, 0, delim='c')
            days = []
            for day in range(len(more_results[0])):
                days.append((day+1)/365.0)
            NPP = np.array(more_results[31])
            cum_NPP = np.cumsum(NPP) * 0.01
            avg_NPP = cum_NPP[-1] / float(days[-1])
            comparison_results[i].append(avg_NPP)
            clean_out(work_path, results_path, ['.png'])
            print
    print
    print "Analyses complete"
    print comparison_results
    c = csv.writer(open(results_path+"productivities.csv", "wb"))
    c.writerow(header)
    for i, comparison in enumerate(comparisons):
        data = [comparison[0], comparison[1], comparison[2]]
        for j in range(len(comparison_results[0])):
            data.append(comparison_results[i][j])
        c.writerow(data)
    print
    print






        # # compute average simulated across sites and plot versus reference
        # length = len(totalc_Mg_ha_set)
        # deadc_Mg_ha_avg = 0
        # totalc_Mg_ha_avg = 0
        # for i in range(length):
        #     if i==0:
        #         deadc_Mg_ha_avg = deadc_Mg_ha_set[0]
        #         totalc_Mg_ha_avg = totalc_Mg_ha_set[0]
        #     else:
        #         deadc_Mg_ha_avg += deadc_Mg_ha_set[i]
        #         totalc_Mg_ha_avg += totalc_Mg_ha_set[i]
        # deadc_Mg_ha_avg /= float(length)
        # totalc_Mg_ha_avg /= float(length)
        #
        #
        #
        # plt.plot(years, deadc_Mg_ha_avg, marker='None', linestyle='-', color='r', label='DayCent simulated dead C')
        # plt.plot(years, totalc_Mg_ha_avg, marker='None', linestyle='-', color='g', label='DayCent simulated total (live+dead) C')
        #
        # ref_years = np.array(reference_values[0])
        # ref_livec = np.array(reference_values[1])
        # ref_deadc = np.array(reference_values[2])
        # ref_totalc = ref_deadc.copy()
        # ref_totalc += ref_livec
        # plt.scatter(ref_years, ref_deadc, marker='o', color='r', label='Greenbook estimated dead C')
        # plt.scatter(ref_years, ref_totalc, marker='o', color='g', label='Greenbook estimated total (live+dead) C')
        #
        # plt.xlabel('Stand age (y)')
        # plt.ylabel('Stand biomass (MgC/ha)')
        # plt.title('%s regional calibration' % region)
        # plt.axhline(0, color='k')
        # plt.legend(loc=2, frameon=False, prop={'size': 11})
        # plt.savefig(results_path+'%s_calibration.png' % region)
        # plt.close()



#########################################################################################

# print the docstring for the user's reference
print
print __doc__

# define key paths
base_path = os.path.dirname(os.path.abspath(__file__))+'/'
sites_path = base_path+"sites/"
work_path = base_path+"workspace/"
library_path = base_path+"input_files/"
ddc_fpath = base_path+"model/DDcentEVI"
ddclist_fpath = base_path+"model/DDClist100"
recent_path = base_path+"current_results/"
archive_path = base_path+"results_archive/"
calibration_path = base_path+"calibration/"

# define model runs in tuple format
#   ( system,   site,        spin-up.bin,                  simulation.sch,              tech,     ax12,  ax5,   ax6,   ax3,   color, linestyle, linewidth, detail, panel8, panel9, panel10, color2, linestyle2, z

# Fig. 2 grouping by FUTURE land use
# run_combos = [
#     ## PASTURE reforestation (unfertilized pine)
#     ('forest', 'NY_coarse', 'pine_pasture_eq_NY.sch',        'forest_NY.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 2, False, False,  True, False, 'saddlebrown', '-', 1),
#     ('forest', 'IA_coarse', 'pine_pasture_eq_IA.sch',        'forest_IA.sch',                '',        True,  False, True,  True, 'saddlebrown', '-', 1, False, False,  True, False, 'saddlebrown', '-', 1),
#     ('forest', 'LA_coarse', 'pine_pasture_eq_LA.sch',        'forest_LA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, False,  True, False, 'saddlebrown', '-', 1),
#
#     ## PASTURE reversion to grassland (unfertilized switchgrass, ungrazed)
#     ('grass',  'NY_coarse', 'grass_pasture_eq.sch',          'grass_NY.sch',              '',        True,  False, False, True, 'g', '-', 2, False, False, True,  False, 'g', '-', 1),
#     ('grass',  'IA_coarse', 'grass_pasture_eq.sch',          'grass_IA.sch',              '',        True,  False, False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),
#     ('grass',  'LA_coarse', 'grass_pasture_eq.sch',          'grass_LA.sch',              '',        True,  False, False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),
#
#     ## PASTURE conversion to energy grass (fertilized switchgrass, ungrazed) - current technology
#     ('bm',     'NY_coarse', 'grass_pasture_eq.sch',          'switchgrass_NY.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 2, False,  True,  False,  False, 'deepskyblue', '-', 2),
#     ('bm',     'IA_coarse', 'grass_pasture_eq.sch',          'switchgrass_IA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True,  False,  False, 'deepskyblue', '-', 2),
#     ('bm',     'LA_coarse', 'grass_pasture_eq.sch',          'switchgrass_LA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True,  False,  False, 'deepskyblue', '-', 2),
#
#     ## PASTURE conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology
#     ('bm_fut', 'NY_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_NY.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 2, False,  True,  False,  False, 'royalblue', '-', 3),
#     ('bm_fut', 'IA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_IA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True,  False,  False, 'royalblue', '-', 3),
#     ('bm_fut', 'LA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_LA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True,  False,  False, 'royalblue', '-', 3),
#
#     ## PASTURE conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology & CCS
#     ('bm_fut', 'NY_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_NY.sch',    'CCS',     True,  False, False, True, 'navy', '-', 2, False,  True,  False,  False, 'navy', '-', 4),
#     ('bm_fut', 'IA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_IA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True,  False,  False, 'navy', '-', 4),
#     ('bm_fut', 'LA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_LA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True,  False,  False, 'navy', '-', 4),
#
#
#     ## CROPLAND reforestation (unfertilized pine)
#     ('forest', 'NY_fine',   'pine_lowSOCcrop_eq_NY.sch',     'forest_NY.sch',                '',        True,  False, True,  True, 'saddlebrown', '-', 2, False, False,  True, False, 'saddlebrown', '-', 1),
#     ('forest', 'IA_fine',   'pine_lowSOCcrop_eq_IA.sch',     'forest_IA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, False,  True, False, 'saddlebrown', '-', 1),
#     ('forest', 'LA_fine',   'pine_lowSOCcrop_eq_LA.sch',     'forest_LA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, False,  True, False, 'saddlebrown', '-', 1),
#
#     ## CROPLAND reversion to grassland (unfertilized switchgrass, ungrazed)
#     ('grass',  'NY_fine',   'grass_lowSOCcrop_eq.sch',       'grass_NY.sch',              '',        True,  False, False, True, 'g', '-', 2, False, False, True,  False, 'g', '-', 1),
#     ('grass',  'IA_fine',   'grass_lowSOCcrop_eq.sch',       'grass_IA.sch',              '',        True,  True,  False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),
#     ('grass',  'LA_fine',   'grass_lowSOCcrop_eq.sch',       'grass_LA.sch',              '',        True,  False, False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),
#
#     ## CROPLAND conversion to energy grass (fertilized switchgrass, ungrazed) - current technology
#     ('bm',     'NY_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_NY.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 2, False,  True,  False,  False, 'deepskyblue', '-', 2),
#     ('bm',     'IA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_IA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True,  False,  False, 'deepskyblue', '-', 2),
#     ('bm',     'LA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_LA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True,  False,  False, 'deepskyblue', '-', 2),
#
#     ## CROPLAND conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology
#     ('bm_fut', 'NY_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_NY.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 2, False,  True,  False,  False, 'royalblue', '-', 3),
#     ('bm_fut', 'IA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_IA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True,  False,  False, 'royalblue', '-', 3),
#     ('bm_fut', 'LA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_LA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True,  False,  False, 'royalblue', '-', 3),
#
#     ## CROPLAND conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology & CCS
#     ('bm_fut', 'NY_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_NY.sch',    'CCS',     True,  False, False, True, 'navy', '-', 2, False,  True,  False,  False, 'navy', '-', 4),
#     ('bm_fut', 'IA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_IA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True,  False,  False, 'navy', '-', 4),
#     ('bm_fut', 'LA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_LA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True,  False,  False, 'navy', '-', 4),
#
#
#     ## FOREST continued growth (unfertilized pine, 70 y.o. stand)
#     ('forest', 'NY_coarse', 'pine_pasture_forest_eq_NY.sch', 'forest_NY.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 2, False, False, False, True,  'saddlebrown', '-', 1),
#     ('forest', 'IA_coarse', 'pine_pasture_forest_eq_IA.sch', 'forest_IA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, False, False, True,  'saddlebrown', '-', 1),
#     ('forest', 'LA_coarse', 'pine_pasture_forest_eq_LA.sch', 'forest_LA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, False, False, True,  'saddlebrown', '-', 1),
#
#     ## FOREST (70 y.o. stand) conversion to energy grass (fertilized switchgrass, ungrazed) - current technology
#     ('forest', 'NY_coarse', 'pine_eq_NY.sch',                'forest_harv_switch_NY.sch', 'current', False, False, False, True, 'deepskyblue', '-',  2, False, False, False, True,  'deepskyblue', '-', 2),
#     ('forest', 'IA_coarse', 'pine_eq_IA.sch',                'forest_harv_switch_IA.sch', 'current', False, False, False, True, 'deepskyblue', '-',  1, False, False, False, True,  'deepskyblue', '-', 2),
#     ('forest', 'LA_coarse', 'pine_eq_LA.sch',                'forest_harv_switch_LA.sch', 'current', False, False, False, True, 'deepskyblue', '-',  1, False, False, False, True,  'deepskyblue', '-', 2),
#
#     ## FOREST (70 y.o. stand) conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology
#     ('forest', 'NY_coarse', 'pine_eq_NY.sch',                'forest_harv_switch_fut_NY.sch', 'mature', False, False, False, True, 'royalblue', '-',  2, False, False, False, True,  'royalblue', '-', 3),
#     ('forest', 'IA_coarse', 'pine_eq_IA.sch',                'forest_harv_switch_fut_IA.sch', 'mature', False, False, False, True, 'royalblue', '-',  1, False, False, False, True,  'royalblue', '-', 3),
#     ('forest', 'LA_coarse', 'pine_eq_LA.sch',                'forest_harv_switch_fut_LA.sch', 'mature', False, False, False, True, 'royalblue', '-',  1, False, False, False, True,  'royalblue', '-', 3),
#
#     ## FOREST (70 y.o. stand) conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology & CCS
#     ('forest', 'NY_coarse', 'pine_eq_NY.sch',                'forest_harv_switch_fut_NY.sch', 'CCS', False, False, False, True, 'navy', '-',  2, False, False, False, True,  'navy', '-', 4),
#     ('forest', 'IA_coarse', 'pine_eq_IA.sch',                'forest_harv_switch_fut_IA.sch', 'CCS', False, False, False, True, 'navy', '-',  1, False, False, False, True,  'navy', '-', 4),
#     ('forest', 'LA_coarse', 'pine_eq_LA.sch',                'forest_harv_switch_fut_LA.sch', 'CCS', False, False, False, True, 'navy', '-',  1, False, False, False, True,  'navy', '-', 4),
#     ]

# Fig. 2 grouping by INITIAL land use, assuming no significant dependency of bioenergy performance on equilibrium land use
run_combos = [
    ## PASTURE reforestation (unfertilized pine)
    ('forest', 'NY_coarse', 'pine_pasture_eq_NY.sch',        'forest_NY.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 2, False, True,  False, False, 'saddlebrown', '-', 1),
    ('forest', 'IA_coarse', 'pine_pasture_eq_IA.sch',        'forest_IA.sch',                '',        True,  False, True,  True, 'saddlebrown', '-', 1, False, True,  False, False, 'saddlebrown', '-', 1),
    ('forest', 'LA_coarse', 'pine_pasture_eq_LA.sch',        'forest_LA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, True,  False, False, 'saddlebrown', '-', 1),

    ## PASTURE reversion to grassland (unfertilized switchgrass, ungrazed)
    ('grass',  'NY_coarse', 'grass_pasture_eq.sch',          'grass_NY.sch',              '',        True,  False, False, True, 'g', '-', 2, False, False, True,  False, 'g', '-', 1),
    ('grass',  'IA_coarse', 'grass_pasture_eq.sch',          'grass_IA.sch',              '',        True,  False, False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),
    ('grass',  'LA_coarse', 'grass_pasture_eq.sch',          'grass_LA.sch',              '',        True,  False, False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),

    ## PASTURE conversion to energy grass (fertilized switchgrass, ungrazed) - current technology
    ('bm',     'NY_coarse', 'grass_pasture_eq.sch',          'switchgrass_NY.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 2, False,  True, True,  False, 'deepskyblue', '-', 2),
    ('bm',     'IA_coarse', 'grass_pasture_eq.sch',          'switchgrass_IA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True, True,  False, 'deepskyblue', '-', 2),
    ('bm',     'LA_coarse', 'grass_pasture_eq.sch',          'switchgrass_LA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True, True,  False, 'deepskyblue', '-', 2),

    ## PASTURE conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology
    ('bm_fut', 'NY_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_NY.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 2, False,  True, True,  False, 'royalblue', '-', 3),
    ('bm_fut', 'IA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_IA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True, True,  False, 'royalblue', '-', 3),
    ('bm_fut', 'LA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_LA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True, True,  False, 'royalblue', '-', 3),

    ## PASTURE conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology & CCS
    ('bm_fut', 'NY_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_NY.sch',    'CCS',     True,  False, False, True, 'navy', '-', 2, False,  True, True,  False, 'navy', '-', 4),
    ('bm_fut', 'IA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_IA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True, True,  False, 'navy', '-', 4),
    ('bm_fut', 'LA_coarse', 'grass_pasture_eq.sch',          'switchgrass_fut_LA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True, True,  False, 'navy', '-', 4),


    ## CROPLAND reforestation (unfertilized pine)
    ('forest', 'NY_fine',   'pine_lowSOCcrop_eq_NY.sch',     'forest_NY.sch',                '',        True,  False, True,  True, 'saddlebrown', '-', 2, False, True,  False, False, 'saddlebrown', '-', 1),
    ('forest', 'IA_fine',   'pine_lowSOCcrop_eq_IA.sch',     'forest_IA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, True,  False, False, 'saddlebrown', '-', 1),
    ('forest', 'LA_fine',   'pine_lowSOCcrop_eq_LA.sch',     'forest_LA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, True,  False, False, 'saddlebrown', '-', 1),

    ## CROPLAND reversion to grassland (unfertilized switchgrass, ungrazed)
    ('grass',  'NY_fine',   'grass_lowSOCcrop_eq.sch',       'grass_NY.sch',              '',        True,  False, False, True, 'g', '-', 2, False, False, True,  False, 'g', '-', 1),
    ('grass',  'IA_fine',   'grass_lowSOCcrop_eq.sch',       'grass_IA.sch',              '',        True,  True,  False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),
    ('grass',  'LA_fine',   'grass_lowSOCcrop_eq.sch',       'grass_LA.sch',              '',        True,  False, False, True, 'g', '-', 1, False, False, True,  False, 'g', '-', 1),

    ## CROPLAND conversion to energy grass (fertilized switchgrass, ungrazed) - current technology
    ('bm',     'NY_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_NY.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 2, False,  True, True,  False, 'deepskyblue', '-', 2),
    ('bm',     'IA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_IA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True, True,  False, 'deepskyblue', '-', 2),
    ('bm',     'LA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_LA.sch',        'current', True,  False, False, True, 'deepskyblue', '-', 1, False,  True, True,  False, 'deepskyblue', '-', 2),

    ## CROPLAND conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology
    ('bm_fut', 'NY_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_NY.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 2, False,  True, True,  False, 'royalblue', '-', 3),
    ('bm_fut', 'IA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_IA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True, True,  False, 'royalblue', '-', 3),
    ('bm_fut', 'LA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_LA.sch',    'mature',  True,  False, False, True, 'royalblue', '-', 1, False,  True, True,  False, 'royalblue', '-', 3),

    ## CROPLAND conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology & CCS
    ('bm_fut', 'NY_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_NY.sch',    'CCS',     True,  False, False, True, 'navy', '-', 2, False,  True, True,  False, 'navy', '-', 4),
    ('bm_fut', 'IA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_IA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True, True,  False, 'navy', '-', 4),
    ('bm_fut', 'LA_fine',   'grass_lowSOCcrop_eq.sch',       'switchgrass_fut_LA.sch',    'CCS',     True,  False, False, True, 'navy', '-', 1, False,  True, True,  False, 'navy', '-', 4),


    ## FOREST continued growth (unfertilized pine, 70 y.o. stand)
    ('forest', 'NY_coarse', 'pine_pasture_forest_eq_NY.sch', 'forest_NY.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 2, False, False, False, True,  'saddlebrown', '-', 1),
    ('forest', 'IA_coarse', 'pine_pasture_forest_eq_IA.sch', 'forest_IA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, False, False, True,  'saddlebrown', '-', 1),
    ('forest', 'LA_coarse', 'pine_pasture_forest_eq_LA.sch', 'forest_LA.sch',                '',        True,  False, False, True, 'saddlebrown', '-', 1, False, False, False, True,  'saddlebrown', '-', 1),

    ## FOREST (70 y.o. stand) conversion to energy grass (fertilized switchgrass, ungrazed) - current technology
    ('forest', 'NY_coarse', 'pine_eq_NY.sch',                'forest_harv_switch_NY.sch', 'current', False, False, False, True, 'deepskyblue', '-',  2, False, False, False, True,  'deepskyblue', '-', 2),
    ('forest', 'IA_coarse', 'pine_eq_IA.sch',                'forest_harv_switch_IA.sch', 'current', False, False, False, True, 'deepskyblue', '-',  1, False, False, False, True,  'deepskyblue', '-', 2),
    ('forest', 'LA_coarse', 'pine_eq_LA.sch',                'forest_harv_switch_LA.sch', 'current', False, False, False, True, 'deepskyblue', '-',  1, False, False, False, True,  'deepskyblue', '-', 2),

    ## FOREST (70 y.o. stand) conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology
    ('forest', 'NY_coarse', 'pine_eq_NY.sch',                'forest_harv_switch_fut_NY.sch', 'mature', False, False, False, True, 'royalblue', '-',  2, False, False, False, True,  'royalblue', '-', 3),
    ('forest', 'IA_coarse', 'pine_eq_IA.sch',                'forest_harv_switch_fut_IA.sch', 'mature', False, False, False, True, 'royalblue', '-',  1, False, False, False, True,  'royalblue', '-', 3),
    ('forest', 'LA_coarse', 'pine_eq_LA.sch',                'forest_harv_switch_fut_LA.sch', 'mature', False, False, False, True, 'royalblue', '-',  1, False, False, False, True,  'royalblue', '-', 3),

    ## FOREST (70 y.o. stand) conversion to energy grass (fertilized switchgrass, ungrazed) - mature technology & CCS
    ('forest', 'NY_coarse', 'pine_eq_NY.sch',                'forest_harv_switch_fut_NY.sch', 'CCS', False, False, False, True, 'navy', '-',  2, False, False, False, True,  'navy', '-', 4),
    ('forest', 'IA_coarse', 'pine_eq_IA.sch',                'forest_harv_switch_fut_IA.sch', 'CCS', False, False, False, True, 'navy', '-',  1, False, False, False, True,  'navy', '-', 4),
    ('forest', 'LA_coarse', 'pine_eq_LA.sch',                'forest_harv_switch_fut_LA.sch', 'CCS', False, False, False, True, 'navy', '-',  1, False, False, False, True,  'navy', '-', 4),
    ]

labels = [
    ('Secondary pine forest growth', 'saddlebrown', '-', 2),
    ('Secondary grassland growth', 'g', '-', 2),
    ('Current switchgrass biofuels', 'deepskyblue', '-', 2),
    ('Mature switchgrass biofuels', 'royalblue', '-', 2),
    ('Mature switchgrass biofuels w/ CCS', 'navy', '-', 2)
    ]

cal_sets = [
    ('NY_WayneCo',
     [[0, 5, 15, 25, 35, 45, 55, 65, 75, 85],
      [2.1, 9.5, 30.4, 46.5, 59.4, 71.1, 80.3, 88.4, 95.9, 102.8],
      [20.4, 16.5, 13.3, 11.4, 10.4, 10.1, 10.1, 10.3, 10.7, 11.1]]),
    ('IA',
     [[0, 5, 15, 25, 35, 45, 55, 65, 75, 85],
      [4.2, 9.3, 18.1, 33.4, 50.7, 66.4, 79.9, 92.3, 101.7, 109.5],
      [17.8, 14.2, 9.9, 9.1, 9.2, 9.8, 10.6, 11.4, 12.1, 12.7]]),
    ('LA',
     [[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85],
      [4.2, 15.5, 27.0, 35.9, 45.5, 55.1, 62.6, 69.5, 76.0, 81.8, 87.5, 92.3, 96.7, 100.7, 104.8, 108.1, 111.5, 114.9],
      [9.2, 8.4, 8.1, 7.8, 7.7, 7.8, 7.9, 8.2, 8.3, 8.5, 8.9, 9.1, 9.3, 9.7, 9.9, 10.0, 10.4, 10.6]])
    ]

# initiate interactive routine
prompt = """
Please select from the following options:
   e = *Execute DayCent simulation ensemble and analyze results
   s = run DayCent model *Spin-up (individually or batch)
   a = *Analyze current simulation results set only
   c = perform pine model *Calibration runs
   p = perform *Productivity comparisons across soil gradient
   q = *Quit
"""
command = ""
while True:
    while True:
        command = raw_input(prompt)
        if command in ('e', 's', 'a', 'c', 'p', 'q'):
            break
        else:
            print "   Please try again!"
    if command == 'e':
        print "   Executing DayCent simulation ensemble..."
        execute(base_path, sites_path, work_path, library_path, recent_path, ddc_fpath,
                ddclist_fpath, run_combos)
        print
        print
        print "   Analyzing current simulation results set..."
        analysis(recent_path, archive_path, run_combos, labels)
        print
        print
    elif command == 's':
        print "   Running DayCent model spin-up..."
        spinup(base_path, sites_path, work_path, library_path, ddc_fpath, ddclist_fpath, run_combos)
        print
        print
    elif command == 'a':
        print "   Analyzing current simulation results set..."
        analysis(recent_path, archive_path, run_combos, labels)
        print
        print
    elif command == 'c':
        print "   Running calibration of pine model..."
        calibration(base_path, cal_sets,  calibration_path, work_path, library_path, ddc_fpath, ddclist_fpath)
        print
        print
    elif command == 'p':
        mukey_fpath = '/data/paustian/AFRI/Ecosystem_dynamics_FROM_NETWORK/soil_exploration/OneidaCoNY_texture_gradient.csv'
        CFARM_db_fpath = '/data/paustian/AFRI/simulations/landscape/Hugoton_34.db'
        print "   Running productivity comparisons across soil gradient in %s..." % mukey_fpath
        soil_comparison(mukey_fpath, CFARM_db_fpath, base_path, archive_path, work_path, library_path, ddc_fpath, ddclist_fpath)
        print
        print
    elif command == 'q':
        print "   Quitting..."
        print
        print
        exit()
