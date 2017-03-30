#!/bin/python

"""Routine to XXX
"""


from db_tools import csv_to_sql
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
import os
import sqlite3


def interactive_soils(soil_atts_fpath):
    """
    """

    def triangule(sand, clay):
        x = (-clay/2.0) - sand  # see 'AFRI grass datasets v011.xlsx'
        y = clay
        return x, y

    def plot_triangle():

        # plot dark triangle outline & labels
        fs = 16
        matplotlib.rcParams.update({'font.size': fs})
        corners = [(0, 0), (-1, 0), (-.5, 1), (0, 0)]
        xy = zip(*corners)
        plt.plot(xy[0], xy[1], color='black', marker=None, linewidth=1.5)
        plt.axis('off')
        plt.text(-0.5, -0.05, ' <-- sand <--', ha='center', va='center', fontsize=fs)
        plt.text(-0.8, 0.5, ' --> clay -->', ha='center', va='center', fontsize=fs, rotation=60)
        plt.text(-0.2, 0.5, ' --> silt -->', ha='center', va='center', fontsize=fs, rotation=-60)

        # plot fine texture class divisions
        soil_class_boundaries = {
            "sand": [(.85, .15, 0), (0.9, 0, 0.1)],
            "loamy sand": [(.7, 0.3, 0), (.85, 0, .15)],
            "silt": [(.2, .8, 0), (.08, .8, .12), (0, .88, .12)],
            "middle": [(.45, 0, .55), (.45, .28, .27), (0, .73, .27)],
            "silt loam": [(.5, .5, 0), (.23, .5, .27)],
            "clay": [(.45, .15, .4), (0, .6, .4)],
            "silts": [(0, .4, .6), (.2, .4, .4), (.2, .53, .27)],
            "sandy clay": [(.65, 0, .35), (.45, .2, .35)],
            "sandy loam": [(.8, 0, .2), (.52, .28, .2), (.45, .28, .27)],
            "loam": [(.52, .28, .2), (.52, .41, .07), (.43, .5, .07)]
        }

        for key in soil_class_boundaries:
            vertex_list = soil_class_boundaries[key]
            vertices_x = []
            vertices_y = []
            for vertex in vertex_list:
                x, y = triangule(vertex[0], vertex[2])
                vertices_x.append(x)
                vertices_y.append(y)
            plt.plot(vertices_x, vertices_y, color='grey', marker=None, linewidth=0.5, zorder=0)

        soil_class_lables = {
            "sand": (.92, .05, .03),
            "loamy\nsand": (.81, .15, .04),
            "sandy\nloam": (.62, .26, .12),
            "sandy\nclay loam": (.6, .12, .28),
            "sandy\nclay": (.52, .07, .41),
            "loam": (.41, .42, .17),
            "clay\nloam": (.33, .33, .34),
            "clay": (.23, .23, .54),
            "silt": (.08, .87, .05),
            "silt\nloam": (.23, .63, .14),
            "silty\nclay": (.06, .47, .47),
            "silty\nclay loam": (.1, .57, .33),
        }

        label_size = 10
        for key in soil_class_lables:
            center = soil_class_lables[key]
            x, y = triangule(center[0], center[2])
            plt.text(x, y, key, color='grey', ha='center', va='center', fontsize=label_size, zorder=0)

    def plot_data(x_set, y_set, facecolors, edgecolor, edgewidth=0.5, point_sizes=16, clickable=False):
        if clickable:
            ax.scatter(x_set, y_set, s=point_sizes, facecolors=facecolors, edgecolors=edgecolor, lw=edgewidth, picker=True)  # 5 points tolerance
        else:
            ax.scatter(x_set, y_set, s=point_sizes, facecolors=facecolors, edgecolors=edgecolor, lw=edgewidth)

    def plot_legend():
        # create dummy entries for legend
        if color_scheme:
            for i, entry in enumerate(legend_increments):
                plt.plot([0, 0], [0, 0], color=legend_colors[i], label=str(legend_increments[i]), linewidth=12)
            plt.legend(bbox_to_anchor=(1.1, 1.1), title=legend_title, prop={'size': 10}, frameon=False)

    def color_brew(value, value_range):
        span = ((value-value_range[0]) / float(value_range[1]-value_range[0]))
        red = 255*(1-span)
        green = 150
        blue = 150
        color = (red/255.0, green/255.0, blue/255.0)
        return color

    def on_pick(event):
        ind = event.ind[0]
        x_pick = xs[ind]
        y_pick = ys[ind]

        # determine the corresponding index that is common across all Xs and Ys, and the associate mukey
        x_indices = [i for i, x in enumerate(xs) if x == x_pick]
        y_indices = [i for i, y in enumerate(ys) if y == y_pick]
        index = list(set(x_indices).intersection(y_indices))[0]
        mukey = mukeys[index]

        # update the soil selection list accordingly, and refresh the plot
        if mukey in selected_mukeys:
            mukey_index = selected_mukeys.index(mukey)
            del selected_mukeys[mukey_index]
            del selected_xs[mukey_index]
            del selected_ys[mukey_index]
            del selected_sizes[mukey_index]
        else:
            selected_mukeys.append(mukey)
            mukey_index = mukeys.index(mukey)
            selected_xs.append(xs[mukey_index])
            selected_ys.append(ys[mukey_index])
            selected_sizes.append(area_weighted_sizes[mukey_index])
        print 'Current soil mukey list:', selected_mukeys
        plt.cla()
        plot_triangle()
        plot_data(xs, ys, colors, 'k', point_sizes=area_weighted_sizes, clickable=True)
        plot_data(selected_xs, selected_ys, 'none', 'r', edgewidth=2, point_sizes=selected_sizes)
        plot_legend()
        plt.show()
        print

    # read soil attributes data to a database format
    print
    print 'Uploading soil attributes data from %s to database...' % soil_atts_fpath
    solution_path, solution_file = os.path.split(soil_atts_fpath)
    db_fpath = solution_path+'/soils.db'
    csv_to_sql(soil_atts_fpath, db_fpath, 'soils', line_block=1000)

    # extract list of unique soils, and associated total areas and n_rate, subject to filtering by land type
    print 'Querying database for unique soils and associated properties...'
    query = "SELECT mukey_run, SUM(hectares), soil_class, avg_sand, avg_silt, avg_clay, depth, total_WHC FROM soils GROUP BY mukey_run"
    print "   Query: %s" % query
    con = sqlite3.connect(db_fpath)
    with con:
        cur = con.cursor()
        cur.execute(query)
        full_soil_list = cur.fetchall()
    os.remove(db_fpath)

    # sort data and transpose for further processing
    print 'Processing soil attributes for plotting as triangle point cloud...'
    print

    filtering = True
    while filtering:
        relative_area_threshold = raw_input("Enter a minimum soil area threshold (in ha) for inclusion in analysis; set to zero to continue without filtering: ")
        filtered_soil_list = []
        total_area = 0
        included_area = 0
        total_soil_counter = 0
        included_soil_counter = 0
        for i, row in enumerate(full_soil_list):
            data = list(row)
            area = data[1]
            total_area += area
            total_soil_counter += 1
            if area > float(relative_area_threshold):
                filtered_soil_list.append(data)
                included_area += area
                included_soil_counter += 1
        print "   You have retained %i/%i soils (%.3f) covering %.0f/%.0f hectares (%.3f). Close plot to continue..." % \
              (included_soil_counter, total_soil_counter, included_soil_counter/float(total_soil_counter),
              included_area, total_area, included_area/float(total_area))
        filtered_soil_list.sort(key=itemgetter(1), reverse=True)   # sort data by area ascending
        mukeys, areas, soil_classes, sands, silts, clays, depths, whcs = zip(*filtered_soil_list)

        # translate textures to triangle plot coordinates
        xs = []
        ys = []
        for i in range(len(mukeys)):
            sand = sands[i]
            clay = clays[i]
            x, y = triangule(sand, clay)
            xs.append(x)
            ys.append(y)

        # normalize & scale areass
        scalar = 12000
        area_weighted_sizes = np.asarray(areas)
        total_area = np.sum(areas)
        area_weighted_sizes /= total_area
        area_weighted_sizes *= scalar

        # plot base data
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plot_triangle()
        plot_data(xs, ys, 'b', 'k', point_sizes=area_weighted_sizes)
        plt.show()
        trigger = raw_input("   Select 'c' to continue, or press enter to filter again: ")
        print
        if trigger == 'c':
            plt.close()
            filtering = False

    # define point colors
    whc_colors = []
    whc_range = (min(whcs), max(whcs))
    whc_span = max(whcs) - min(whcs)
    for i, whc in enumerate(whcs):
        whc_color = color_brew(whc, whc_range)
        whc_colors.append(whc_color)

    depth_colors = []
    depth_range = (min(depths), max(depths))
    depth_span = max(depths) - min(depths)
    for i, depth in enumerate(depths):
        depth_color = color_brew(depth, depth_range)
        depth_colors.append(depth_color)

    # prompt user for inputs specifying point size and color display
    colors = 'b'
    legend_increments = []
    legend_colors = []
    legend_title = ''
    color_scheme = raw_input("For point coloring by depth type 'd'; by WCH type 'w', otherwise press enter to continue with unifor point color: ")
    print
    if color_scheme == 'd':
        colors = depth_colors
        legend_increments = [min(depths), min(depths)+0.25*depth_span, min(depths)+0.5*depth_span, min(depths)+0.75*depth_span, max(depths)]
        legend_colors = [color_brew(n, depth_range) for n in legend_increments]
        legend_title = 'Total soil column depth (cm)'
    elif color_scheme == 'w':
        colors = whc_colors
        legend_increments = [min(whcs), min(whcs)+0.25*whc_span, min(whcs)+0.5*whc_span, min(whcs)+0.75*whc_span, max(whcs)]
        legend_colors = [color_brew(n, whc_range) for n in legend_increments]
        legend_title = 'Total soil column WHC (cm)'

    # create data structure for storing selected soils and associated attributes, replot soils, and enter interactive mode
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_triangle()
    plot_data(xs, ys, colors, 'k', point_sizes=area_weighted_sizes, clickable=True)
    plot_legend()
    selected_mukeys = []
    selected_xs = []
    selected_ys = []
    selected_sizes = []
    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()


interactive_soils('/Users/johnfield/Desktop/local_python_code/Interactive_plots/36117_atts_soils.csv')