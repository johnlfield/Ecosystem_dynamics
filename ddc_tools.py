#!/bin/python

"""Collection of functions for analyzing manipulating DayCent input files. 
"""

import csv
import glob
import math
import os


def fcs_wps(mukey, layers=''):
    #read appropriate soils.in content to a python list
    mukey = str(mukey)
    soil_path = "/data/paustian/ernie/SSURGO_master_script/soil_test2/"
    soil_fpath = soil_path+mukey[:-3]+"/"+mukey+".in"
    cont = [[]]
    data_input = open(soil_fpath, 'r')
    for line in data_input:
        cont.append(line.split())
    del cont[0]

    #convert all entries in the 2D list to float format where possible, or zero in the case
    #of very small numbers recorded in scientific notation
    for k in range(len(cont)):
        for l in range(len(cont[k])):
            cont[k][l] = float(cont[k][l])

    # extract all field capacity and wilting point values
    fcs = []
    wps = []
    for i in range(len(cont)):
        if not layers:
            fc = cont[i][3]
            wp = cont[i][4]
            fcs.append(fc)
            wps.append(wp)
        else:
            if 1+i <= layers:
                fc = cont[i][3]
                wp = cont[i][4]
                fcs.append(fc)
                wps.append(wp)

    return fcs, wps


def whc_tot(mukey, layers=''):
    """Reads the DayCent-formatted soils.in file for the specified SSURGO mukey, and
    computes the corresponding total water holding capacity for the full soil profile
    in cm. 
    Args:
        mukey- soil SSURGO map unit key (int)
    """
    #read appropriate soils.in content to a python list
    mukey = str(mukey)
    soil_path = "/data/paustian/ernie/SSURGO_master_script/soil_test2/"
    soil_fpath = soil_path+mukey[:-3]+"/"+mukey+".in"
    cont = [[]]
    data_input = open(soil_fpath, 'r')
    for line in data_input:
        cont.append(line.split())
    del cont[0]

    #convert all entries in the 2D list to float format where possible, or zero in the case
    #of very small numbers recorded in scientific notation
    for k in range(len(cont)):
        for l in range(len(cont[k])):
            cont[k][l] = float(cont[k][l])

    #loop through list and compute the water holding capacity increment represented in 
    #each line
    min_h2o_evap = 0
    min_h2o = 0
    max_h2o = 0
    whc = 0
    for i in range(len(cont)):
        if not layers:
            depth = cont[i][1] - cont[i][0]
            FC = cont[i][3]
            WP = cont[i][4]
            WHC = FC - WP
            if i != 0:
                min_h2o_evap += depth*WP
            min_h2o += depth*WP
            max_h2o += depth*FC
            whc += depth*WHC
        else:
            if 1+i <= layers:
                depth = cont[i][1] - cont[i][0]
                FC = cont[i][3]
                WP = cont[i][4]
                WHC = FC - WP
                if i != 0:
                    min_h2o_evap += depth*WP
                min_h2o += depth*WP
                max_h2o += depth*FC
                whc += depth*WHC
    if layers:
        if layers > len(cont):
            print "NOTE: specified layer limit exceeds number of layers found in soils.in file"

    return whc, min_h2o, max_h2o


def avg_text(mukey, layers):
    """Reads the DayCent-formatted soils.in file for the specified SSURGO mukey, and
    computes the corresponding depth-weighted average soil texture.
    Args:
        mukey- soil SSURGO map unit key (int)
    """
    #read appropriate soils.in content to a python list
    mukey = str(mukey)
    soil_path = "/data/paustian/ernie/SSURGO_master_script/soil_test2/"
    soil_fpath = soil_path+mukey[:-3]+"/"+mukey+".in"
    cont = [[]]
    data_input = open(soil_fpath, 'r')
    for line in data_input:
        cont.append(line.split())
    del cont[0]

    #convert all entries in the 2D list to float format where possible, or zero in the case
    #of very small numbers recorded in scientific notation
    for k in range(len(cont)):
        for l in range(len(cont[k])):
            cont[k][l] = float(cont[k][l])

    #loop through list and compute the depth-weighted fraction of each texture component
    sand_tot = 0
    silt_tot = 0
    clay_tot = 0
    for i in range(len(cont)):
        if i+1 <= layers:
            depth = float(cont[i][1]) - float(cont[i][0])
            sand = float(cont[i][7])
            clay = float(cont[i][8])
            silt = round(1-sand-clay, 2)
            sand_tot += sand * depth
            silt_tot += silt * depth
            clay_tot += clay * depth
            final_depth = float(cont[i][1])

    if layers > len(cont):
        print "NOTE: specified layer limit exceeds number of layers found in soils.in file"

    # normalize by total depth
    sand_avg = sand_tot/final_depth
    silt_avg = silt_tot/final_depth
    clay_avg = clay_tot/final_depth

    return sand_avg, silt_avg, clay_avg, final_depth


def surf_text(mukey):
    """Reads the DayCent-formatted soils.in file for the specified SSURGO mukey, and
    returns the sand/silt/clay fractions for the top layer. 
    Args:
        mukey- soil SSURGO map unit key (int)
    """
    #read appropriate soils.in content to a python list
    mukey = str(mukey)
    soil_path = "/data/paustian/ernie/SSURGO_master_script/soil_test2/"
    soil_fpath = soil_path+mukey[:-3]+"/"+mukey+".in"
    data_input = open(soil_fpath, 'r')
    top = next(data_input).split()
    sand = float(top[7])
    clay = float(top[8])
    silt = round(1-sand-clay, 2)

    return sand, silt, clay


def replace_100(in_fpath, replacements, out_fpath):
    """Generic routine to read a template .100 file and overwrite wildcards with values.
    :param in_fpath:  path to template file
    :param replacements:  dictionary of replacement values (strings)
    :param out_fpath:  output file path
    :return:
    """
    reader = open(in_fpath)
    writer = open(out_fpath, 'w')
    for row in reader:
        for src, target in replacements.iteritems():
            row = row.replace(src, target)
        writer.write(row)
    reader.close()
    writer.close()


def param_lookup(param_fpath, param_name, version=-999):
    """Fill out the specified .100 file templates based on parameters read from the specified parameter .csv file.
    :param param_fpath:  .csv-format run parameter file
    :param template_path:  path containing all *.100.template files with wildcard strings
    :param target_fpath:  path where resulting *.100 files will be saved
    :return:
    """

    # save the column names and read the last line of parameters in a params.csv file
    params_csv_object = open(param_fpath, 'rU')
    lines = csv.reader(params_csv_object)
    keys = lines.next()   # save column names as keys
    for line in lines:
        if line[3]:   # check that data is present, any extra blank lines saved at end of file won't override real data
            params = line
            read_version = int(params[0])
        # break out of the loop once the target parameter version have been read, if specified
        if read_version == version:
            break
    params_csv_object.close()
    print "Updating .100 files based on parameter version %i..." % read_version

    # propagate parameter data into dictionary
    dictionary = {}
    for i, key in enumerate(keys):
        dictionary[key] = params[i]

    return dictionary[param_name]


def param_100(param_fpath, template_path, target_path, version=-999):
    """Fill out the specified .100 file templates based on parameters read from the specified parameter .csv file.
    :param param_fpath:  .csv-format run parameter file
    :param template_path:  path containing all *.100.template files with wildcard strings
    :param target_fpath:  path where resulting *.100 files will be saved
    :return:
    """

    # save the column names and read the last line of parameters in a params.csv file
    params_csv_object = open(param_fpath, 'rU')
    lines = csv.reader(params_csv_object)
    keys = lines.next()   # save column names as keys
    for line in lines:
        if line[3]:   # check that data is present, any extra blank lines saved at end of file won't override real data
            params = line
            read_version = int(params[0])
        # break out of the loop once the target parameter version have been read, if specified
        if read_version == version:
            break
    params_csv_object.close()
    print "Updating .100 files based on parameter version %i..." % read_version

    # propagate parameter data into dictionary and perform replacement
    dictionary = {}
    for i, key in enumerate(keys):
        dictionary[key] = params[i]

    # iterate through and perform dictionary-based replacement for each .template file found
    for fpath in glob.glob(os.path.join(template_path, "*")):
        if fpath.endswith(".template"):
            template_fpath = fpath
            target_file = fpath.split("/")[-1][:-9]
            target_fpath = target_path+target_file
            print "   Generating parameterized %s" % target_fpath
            replace_100(template_fpath, dictionary, target_fpath)

    return read_version


def soil_h2o_detail(vswc_fpath, mukey, layer_lim, cwscoef11, cwscoef12, start_year, end_year):
    """Routine for summarizing visualizing seasonal soil moisture dynamics from a vswc.out file.  The routine sums up
    daily total root zone volumetric water content, and plots it relative to water holding capacity limits for the
    wettest and driest year in the file.
    :param vswc_fpath:  full path to vswc.out for analysis
    :param mukey:  associated soils.in map unit key for determination of
    :param layers:  number of layers eligible for root development, should be consistent with crop.100 CLAYPG
    :return:  none
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # define soils.in layer depth structure and layer water holding capacities
    depths = [2, 3, 5, 10, 10, 15, 15, 15, 15, 15, 15, 30, 30, 30]
    fcs, wps = fcs_wps(mukey, layers=layer_lim)

    # import file contents
    vswc_file = open(vswc_fpath, 'rU')
    contents = []
    for layer in vswc_file:
        entries = layer.split()
        for i in range(len(entries)):
            entries[i] = float(entries[i])

        # convert the year from float to int, and note the day
        year = int(math.floor(entries[0]))
        day = int(entries[1])

        # sum depth-weighted moisture content for layers within limit
        water_vol = 0
        for j in range(len(entries)):
            if j+1 <= layer_lim:
                water_vol += depths[j] * entries[j+2]

        # determine layer with maximum relative water content, and associated value
        rwcs = []  # relative water content
        for j in range(layer_lim):
            vwc = entries[j+2]
            rwc = (vwc-wps[j]) / (fcs[j]-wps[j])
            rwcs.append(rwc)
        max_rwc = max(rwcs)

        contents.append([year, day, water_vol, max_rwc])

    # # identify wettest and driest year w/in permissible range
    # annual_totals = [0]
    # years = []
    # latest_year = 200000
    # j = 0
    # for i in range(len(contents)):
    #     year = int(math.floor(contents[i][0]))
    #     if year <= latest_year:
    #         annual_totals[j] += contents[i][2]
    #     else:
    #         annual_totals.append(0)
    #         years.append(contents[i][0])
    #         j += 1
    #     latest_year = year
    # year2 = years[0]
    # years.insert(0, year2-1)
    # perm_years = []
    # perm_annual_totals = []
    # for year in years:
    #     if start_year <= year <= end_year:
    #         perm_years.append(year)
    #         perm_annual_totals.append(annual_totals[years.index(year)])
    # wet_year = years[annual_totals.index(max(annual_totals))]
    # dry_year = years[annual_totals.index(min(annual_totals))]
    wet_year = start_year
    dry_year = end_year

    # extract wettest and driest year series within possible range
    wet_svwc = []
    dry_svwc = []
    wet_days = []
    dry_days = []
    wet_rwcs = []
    dry_rwcs = []
    for i in range(len(contents)):
        year = int(math.floor(contents[i][0]))
        if year == wet_year:
            wet_days.append(contents[i][1])
            wet_svwc.append(contents[i][2])
            wet_rwcs.append(contents[i][3])
        if year == dry_year:
            dry_days.append(contents[i][1])
            dry_svwc.append(contents[i][2])
            dry_rwcs.append(contents[i][3])

    # plot
    work_path = os.path.dirname(vswc_fpath)
    scenario = vswc_fpath.split('/')[-1][:-4]
    png_fpath = work_path+'/%s.png' % scenario
    f, axes = plt.subplots(2, 2, sharex=True)
    axes[0][0].plot(wet_days, wet_svwc, zorder=2)
    axes[0][0].set_xlim(100, 300)
    axes[0][0].set_ylim(0, 30)
    axes[0][0].set_title('%i' % wet_year)
    axes[0][0].set_ylabel('Root zone\nwater content (cm)')
    axes[0][1].plot(dry_days, dry_svwc, zorder=2)
    axes[0][1].set_xlim(100, 300)
    axes[0][1].set_ylim(0, 30)
    axes[0][1].set_title('%i' % dry_year)

    axes[1][0].plot(wet_days, wet_rwcs, zorder=2)
    axes[1][0].fill_between(wet_days, 0, wet_rwcs, color='c', zorder=0)
    axes[1][0].set_ylabel('Maximum layer\nrelative water content')
    axes[1][1].plot(dry_days, dry_rwcs, zorder=2)
    axes[1][1].fill_between(dry_days, 0, dry_rwcs, color='c', zorder=0)
    axes[1][0].set_ylim(0, 1)
    axes[1][1].set_ylim(0, 1)

    # add limits & save
    whc, min_h2o, max_h2o = whc_tot(mukey, layers=layer_lim)
    production_levels = [0.25, 0.5, 0.75]
    colors = ['red', 'orange', 'green']
    for i in range(len(production_levels)):
        rwcf = cwscoef11 - (math.log((1-production_levels[i])/production_levels[i]) / cwscoef12)
        axes[1][0].axhline(y=rwcf, color=colors[i], zorder=1)
        axes[1][1].axhline(y=rwcf, color=colors[i], zorder=1)
        axes[1][1].plot(200, 0, ms='none', color=colors[i], label='%s NPP' % str(production_levels[i]))
    axes[1][1].legend(loc='upper right', prop={'size': 10})

    axes[0][0].axhline(y=min_h2o, color='red', zorder=1)
    axes[0][0].axhline(y=max_h2o, color='green', zorder=1)
    axes[0][1].axhline(y=min_h2o, color='red', zorder=1)
    axes[0][1].axhline(y=max_h2o, color='green', zorder=1)
    axes[0][1].plot(200, 0, ms='none', color='green', label='field capacity')
    axes[0][1].plot(200, 0, ms='none', color='red', label='wilting point')
    axes[0][1].legend(loc='upper right', prop={'size': 10})

    f.text(0.5, 0.02, 'Day of year', ha='center')
    plt.savefig(png_fpath)
    plt.close()


def soil_class(sand, silt, clay):

    texture = 'ERROR'
    if sand+silt+clay != 1:
        texture = 'DOES NOT SUM TO UNITY'
    elif (silt*100 + 1.5*clay*100) < 15:
        texture = "sand"
    elif (silt*100 + 1.5*clay*100) >= 15 and (silt*100 + 2*clay*100) < 30:
        texture = "loamy sand"
    elif ((silt*100+2*clay*100) >= 30 and clay*100 >= 7 and clay*100 < 20 and sand*100 > 52) or ((silt*100+2*clay*100) >= 30 and silt*100 < 50 and clay*100 < 7):
        texture = "sandy loam"
    elif sand*100 <= 52 and silt*100 >= 28 and silt*100 < 50 and clay*100 >= 7 and clay*100 < 27:
        texture = "loam"
    elif (clay*100 >= 12 and silt*100 >= 50 and clay*100 < 27) or (silt*100 >= 50 and silt*100 < 80 and clay*100 < 12):
        texture = "silt loam"
    elif silt*100 >= 80 and clay*100 < 12:
        texture = "silt"
    elif sand*100 > 45 and silt*100 < 28 and clay*100 >= 20 and clay*100 < 35:
        texture = "sandy clay loam"
    elif sand*100 <= 20 and clay*100 >= 27 and clay*100 < 40:
        texture = "silty clay loam"
    elif sand*100 > 20 and sand*100 <= 45 and clay*100 >= 27 and clay*100 < 40:
        texture = "clay loam"
    elif sand*100 > 45 and clay*100 >= 35:
        texture = "sandy clay"
    elif silt*100 >= 40 and clay*100 >= 40:
        texture = "silty clay"
    elif sand*100 <= 45 and silt*100 < 40 and clay*100 >= 40:
        texture = "clay"
    return texture