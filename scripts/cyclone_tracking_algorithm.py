#!/usr/bin/env python
# coding: utf-8

# =============================================================================

# Cyclone tracking algorithm, as used in Chen et al. (2025), adjusted for
#   compatibility with CRCM5 rotated coordinate system

# =============================================================================

# Associated publication: Weaver, M.M. et al. (2026). Historical nor'easter
#                         intensity trends in reanalysis and regional climate
#                         simulations: A cross-dataset perspective.

# =============================================================================

# This file includes all functions written by Chen et al. (2025), some of which
#   are not used in the present work. New and/or edited functions were added
#   to Chen et al. (2025)'s algorithm to make the code compatible with non-ERA5
#   data structures.

# =============================================================================

# Import necessary packages

import time, os, sys, scipy
import numpy as np
import cartopy as ct
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')  # Use non-interactive backend suitable for headless environments
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.animation as animation
import skimage as ski
from datetime import datetime, timedelta
import cftime
from pyproj import CRS, Transformer
from shapely.geometry import Point, Polygon

# =============================================================================

# Plotting utility functions (not used in the present work)

# plots filled contours
def plot_function(data_arr, lons, lats, clevels):
    ax = plt.axes(projection=ct.crs.PlateCarree())
    ax.set_extent([-90, right_lon, 20, top_lat], crs=ct.crs.PlateCarree())
    ax.coastlines()
    plt.contourf(lons, lats, data_arr, levels=clevels, cmap='jet', 
                 transform=ct.crs.PlateCarree())
    plt.contour(lons, lats, data_arr, levels=clevels, 
                transform=ct.crs.PlateCarree(), linewidths=0.75)
    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(clevels, cmap.N)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, 
                        orientation='vertical', 
                        label="Mean sea level pressure [Pa]", 
                        fraction=0.031, pad=0.04)
    return ax


# plots filled contours + 1 set of points
def plot_function_2(data_arr, lons, lats, clevels, xcs, ycs):
    ax = plt.axes(projection=ct.crs.PlateCarree())
    ax.set_extent([-90, right_lon, 20, top_lat], crs=ct.crs.PlateCarree())
    ax.coastlines()
    plt.contourf(lons, lats, data_arr, levels=clevels, cmap='jet', 
                 transform=ct.crs.PlateCarree())
    plt.contour(lons, lats, data_arr, levels=clevels, 
                transform=ct.crs.PlateCarree(), linewidths=0.75)
    lonsc1, latsc1 = lat_lon(xcs, ycs, lons, lats)
    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(clevels, cmap.N)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                        ax=ax, orientation='vertical', 
                        label="Mean sea level pressure [Pa]", 
                        fraction=0.031, pad=0.04)
    ax.scatter(lonsc1, latsc1, s=60, facecolors='r', edgecolors='white')
    plt.show()
    return ax


# plots filled contours + 2 sets of points
def plot_function_3(data_arr, lons, lats, clevels, xcs, ycs, xps, yps):
    ax = plt.axes(projection=ct.crs.PlateCarree())
    ax.set_extent([-90, right_lon, 20, top_lat], crs=ct.crs.PlateCarree())
    ax.coastlines()
    plt.contourf(lons, lats, data_arr, levels=clevels, cmap='jet', 
                 transform=ct.crs.PlateCarree())
    plt.contour(lons, lats, data_arr, levels=clevels, 
                transform=ct.crs.PlateCarree(), linewidths=0.75)
    lonsc1, latsc1 = lat_lon(xcs, ycs, lons, lats)
    lonsp1, latsp1 = lat_lon(xps, yps, lons, lats)
    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(clevels, cmap.N)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), 
                        ax=ax, orientation='vertical', 
                        label="Mean sea level pressure [Pa]", 
                        fraction=0.031, pad=0.04)
    ax.scatter(lonsc1, latsc1, s=60, facecolors='r', alpha=0.5)
    ax.scatter(lonsp1, latsp1, s=60, facecolors='b', alpha=0.5)
    plt.show()
    return ax


# plots all storms from dataframe in 6-hour timesteps
def plot_tracks(df_out, datai, timesi, lons, lats, contours):
    xcs_arr, ycs_arr, tarri = get_storm_coords(df_out, timesi, 'c')
    xps_arr, yps_arr, tarri = get_storm_coords(df_out, timesi, 'p')
    for j in range(len(tarri)):
        print('Storm ', j)
        for i in range(len(tarri[j])-1):
            plot_function_3(datai[tarri[j][i]], lons, lats, contours, 
                            xcs_arr[j][i], ycs_arr[j][i], xps_arr[j][i], 
                            yps_arr[j][i])
        print('')


# plots objects from get_objects_mslp output
def plot_obj(t, data, times, lons, lats, contours): # removed lati, loni from Chen et al. (2025)
    df_p = get_objects_mslp(t, data, lons, lats, times, True) # removed lati, loni
    plot_function_3(data[t], lons, lats, contours, df_p['cx_c'].to_numpy(), 
                    df_p['cy_c'].to_numpy(), df_p['cx_p'].to_numpy(), 
                    df_p['cy_p'].to_numpy())
    return df_p


# make .gif for a storm; MSLP map with centers
def make_gif(data, lons, lats, clevels, xcs, ycs, xps, yps, tarri, idi):
    num_frames = len(tarri) - 1
    frame_list = range(num_frames)
    fig = plt.figure()
    ax = plt.axes(projection=ct.crs.PlateCarree())
    ax.set_extent([-90, right_lon, 20, top_lat], crs=ct.crs.PlateCarree())
    ax.coastlines()
    cs = plt.contourf(lons, lats, data[tarri[0]], levels=clevels, cmap='jet', 
                      transform=ct.crs.PlateCarree())
    lonsc1, latsc1 = lat_lon(xcs[0], ycs[0], lons, lats)
    lonsp1, latsp1 = lat_lon(xps[0], yps[0], lons, lats)
    cmap = mpl.cm.jet
    norm = mpl.colors.BoundaryNorm(clevels, cmap.N)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, 
                        orientation='vertical', 
                        label="Mean sea level pressure [Pa]", 
                        fraction=0.031, pad=0.04)
    ax.scatter(lonsc1, latsc1, s=60, facecolors='r', alpha=0.5)
    ax.scatter(lonsp1, latsp1, s=60, facecolors='b', alpha=0.5)

    def animate(frame):
        cs = plt.contourf(lons, lats, data[tarri[frame]], levels=clevels, 
                          cmap='jet', transform=ct.crs.PlateCarree())
        lonsc1f, latsc1f = lat_lon(xcs[frame], ycs[frame], lons, lats)
        lonsp1f, latsp1f = lat_lon(xps[frame], yps[frame], lons, lats)
        ax.scatter(lonsc1f, latsc1f, s=60, facecolors='r', alpha=0.5)
        ax.scatter(lonsp1f, latsp1f, s=60, facecolors='b', alpha=0.5)
        #print(f'Frame: {frame+1}/{num_frames}')
        plt.title(f'{idi}, + {frame*dt} hr')
        return cs

    ofile = idi + '.gif'
    anim = animation.FuncAnimation(fig, animate, frames=frame_list)
    anim.save(ofile)
    plt.close()


# make .gif for all storms in dataframe
def make_gif_all(df_out, data, times, lons, lats, clevels):
    num_storms = len(df_out)
    xcs_arr, ycs_arr, tarri = get_storm_coords(df_out, times, 'c')
    xps_arr, yps_arr, tarri = get_storm_coords(df_out, times, 'p')
    id_list = df_out['ID'].to_numpy()
    for i in range(num_storms):
        make_gif(data, lons, lats, clevels, xcs_arr[i], ycs_arr[i], xps_arr[i], 
                 yps_arr[i], tarri[i], id_list[i])

# =============================================================================

# Storm tracking utility functions

# given a gridpt (xval, yval), returns a box around the gridpt
def grow_box(xval, yval, x_rad, y_rad, data_arr): #, lati, loni):
    x1 = xval-x_rad
    x2 = xval+x_rad
    y1 = yval-y_rad
    y2 = yval+y_rad
    if x1<0:
        x1 = 0
    if x2>len(data_arr[0]):
        x2 = len(data_arr[0])
    if y1<0:
        y1 = 0
    if y2>len(data_arr):
        y2 = len(data_arr)
    borders = [x1, x2, y1, y2]
#    a1 = x2+1
#    a2 = y2+1
#    if x2 >= loni:
#        a1 = loni+1
#    elif y2 >= lati:
#        a2 = lati+1
    return np.arange(x1, x2), np.arange(y1, y2), borders


# converts gridpts to lat/lon coordinates
def lat_lon(gridx, gridy, lons, lats):
    lon = lons[gridy, gridx]
    if lon > 180:
        lon -= 360
    lat = lats[gridy, gridx]
    return lon, lat

# given a list of gridpts, convert to lat/lon coordinates
def grids_to_map(xs, ys, lons, lats):
    lons_coords = [lat_lon(x, 0, lons, lats)[0] for x in xs]
    lats_coords = [lat_lon(0, y, lons, lats)[1] for y in ys]
    return lons_coords, lats_coords


# for a given latitude, converts physical distance to a distance in gridpts
def grid_rad_lon(lat, dist):
    dist_25 = (np.pi/180)*6378*np.cos(np.radians(lat))/4
    x_grid = dist/dist_25
    return np.round(abs(x_grid))


# calculates physical distance between lat/lon coordinates
def phys_dist(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295
    a = 0.5 - np.cos((lat2 - lat1)*p)/2 + np.cos(lat1*p)*np.cos(lat2*p)*(1 - np.cos((lon2 - lon1)*p))/2
    return 12742 * np.arcsin(np.sqrt(a))


# returns list of values from initialized storm track dictionary
def key_val(df_obj):
    dlist = []
    for i in range(len(df_obj)):
        val = {}
        val['cx_p'] = [df_obj['cx_p'][i]]
        val['cy_p'] = [df_obj['cy_p'][i]]
        val['c_p (Pa)'] = [df_obj['c_p (Pa)'][i]]
        #val['x_contour'] = [df_obj['x_contour'][i]]
        #val['y_contour'] = [df_obj['y_contour'][i]]
        val['cx_c'] = [df_obj['cx_c'][i]]
        val['cy_c'] = [df_obj['cy_c'][i]]
        dlist.append(val)
    return dlist


# given storm ID, returns the start time, start time index, & array of time indices
def t_from_id(id, coord_len, timesi, dt):
    mdy = id.split('-')[0]
    t0 = datetime(year=int(mdy[4:]), day=int(mdy[2:4]), 
                  month=int(mdy[:2]), hour=int(id.split('-')[1]))
    t0i = int((t0 - timesi[0]).total_seconds()/3600)
    ts = np.arange(t0i, t0i + dt*coord_len, 1)
    return t0, t0i, ts


# returns gridcoords for a dictionary of storm tracks
def get_storm_coords(df_out, timesi, pc):
    num_storms = len(df_out)
    ids = df_out['ID'].to_numpy()
    if pc == 'c':
        xs = df_out['cx_c'].to_numpy()
        ys = df_out['cy_c'].to_numpy()
    elif pc == 'p':
        xs = df_out['cx_p'].to_numpy()
        ys = df_out['cy_p'].to_numpy()
    xs_arr = []
    ys_arr = []
    tarri = []
    for i in range(num_storms):
        #t0, t0i, ts = t_from_id(ids[i], len(xs[i]), timesi)
        #tarri.append(ts)
        for j in xs[i][:-1]:
            xs_arr.append(int(j))
        for k in ys[i][:-1]:
            ys_arr.append(int(k))
    return xs_arr, ys_arr#, tarri


# prints out storm track dictionary (coords)
def print_tracks(track_d, pc):
    for k, v in track_d.items():
        print(k)
        if pc == 'p':
            print('center_x_pressure: ', v['cx_p'])
            print('center_y_pressure: ', v['cy_p'])
        elif pc == 'c':
            print('center_x_contour: ', v['cx_c'])
            print('center_y_contour: ', v['cy_c'])
        print()


# Interpolation function to define the boundary lines
def interpolate_boundary(start_lat, start_lon, end_lat, end_lon, num_points=100):
    # Linear interpolation between two points
    lats = np.linspace(start_lat, end_lat, num_points)
    lons = np.linspace(start_lon, end_lon, num_points)
    return lats, lons

# Function to check if a point is inside a polygon using Shapely
def is_point_in_polygon(lat_point, lon_point, lat_boundary, lon_boundary):
    polygon = Polygon(zip(lon_boundary, lat_boundary))
    point = Point(lon_point, lat_point)
    return polygon.contains(point)

# determines if storm track passes through the US east coast lat/lon box
# returns the number of times that storm passes within box
def pass_USEC(xps, yps, xcs, ycs, lats, lons):
    passed = 0
    #storm_lon_lat = []
    for i in range(len(xcs)):
        # Handle 'f' as end of storm and skip this track point
        if xcs[i] == 'f' or ycs[i] == 'f' or xps[i] == 'f' or yps[i] == 'f':
            print(f"Encountered 'f' at index {i}, breaking out of loop.")
            break
        # Convert grid indices to lat/lon for each point
        lonp1, latp1 = lat_lon(xps[i], yps[i], lons, lats)
        lonc1, latc1 = lat_lon(xcs[i], ycs[i], lons, lats)
        # Checkpoint: coordinates
        print('Storm lon: ' + str(lonp1))
        print('Storm lat: ' + str(latp1))

        # Check if the point (either p or c) is inside the box (using lat/lon bounds)
        p_in = (usne_W <= lonp1 <= usne_E) and (usne_S <= latp1 <= usne_N)
        c_in = (usne_W <= lonc1 <= usne_E) and (usne_S <= latc1 <= usne_N)
        print(str(p_in))
        print(str(c_in))

        # If either point is inside the box, count it
        if p_in or c_in:
            passed += 1

    return passed


def plot_storm_track(track_coords, storm_id):
    lons, lats = zip(*track_coords)

    plt.figure(figsize=(8, 6))
    plt.plot(lons, lats, marker='o', label=f'Storm {storm_id}')
    # Draw the East Coast box
    plt.plot([usne_W, usne_E, usne_E, usne_W, usne_W],
             [usne_S, usne_S, usne_N, usne_N, usne_S], 'r--', label='USEC Box')

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Storm {storm_id} Track with USEC Box')
    plt.grid()
    plt.legend()

    os.makedirs("storm_plots", exist_ok=True)
    filepath = f"storm_plots/storm_{storm_id}_track.png"
    plt.savefig(filepath)
    plt.close()
    print(f"Saved plot for storm {storm_id} to {filepath}")


# identifies storms that pass through a bounding box around the US east coast
def filter_USEC(df_track, lats, lons):
    passed = []
    num_storms = len(df_track)
    for i in range(num_storms):
        xps = df_track['cx_p'].to_numpy()[i]
        yps = df_track['cy_p'].to_numpy()[i]
        xcs = df_track['cx_c'].to_numpy()[i]
        ycs = df_track['cy_c'].to_numpy()[i]
        storm_passed = pass_USEC(xps, yps, xcs, ycs, lats, lons)
        passed.append(storm_passed)

        # Plot every storm
        #plot_storm_track_with_raw(xps, yps, lons, lats, converted_track, i)
    df_track['passed_USEC'] = passed
    return df_track


# takes storm tracking dataframe and drops storms that don't reach a min SLP
def SLP_filter(df_out, minSLP):
    df_out = df_out.drop(df_out[df_out['min_p (Pa)'] > minSLP].index)
    df_out = df_out.reset_index(drop=True)
    return df_out


# takes storm tracking dataframe and drops storms that don't last 24 hrs
def filter_time(df_out):
    df_out = df_out.drop(df_out[df_out['lifetime (h)'] < 24].index)
    df_out = df_out.reset_index(drop=True)
    return df_out


# takes storm tracking dataframe and drops storms that don't travel
# at least 1000 km (pressure and spatial center)
def filter_dist(df_out):
    df_out = df_out.drop(df_out[df_out['tot_dist_p (km)'] < 1000].index)
    df_out = df_out.drop(df_out[df_out['tot_dist_c (km)'] < 1000].index)
    df_out = df_out.reset_index(drop=True)
    return df_out


# takes storm tracking dataframe and drops storms that form too far west and/or
#     too far southmeant to filter out TCs (south) and Miller B nor'easters 
#     (west), if desired
def filter_gen(df_out, max_longitude, max_latitude, lats, lons):
    # Check the first position of each storm and drop those that don't meet the 
    # longitude and latitude criteria
    mask = df_out.apply(
        lambda row: lat_lon(row['cx_p'][0], row['cy_p'][0], lons, 
                            lats)[0] <= max_longitude or  # Changed < to <=
                    lat_lon(row['cx_p'][0], row['cy_p'][0], lons, 
                            lats)[1] < max_latitude, axis=1
    )
    # Drop the rows where the mask is True
    df_out = df_out[~mask].reset_index(drop=True)
    return df_out

# given storm track dictionary, does some additional analysis,
# then returns a dataframe and writes to .csv
# includes NaNs debugging (distance)
def storm_track_out(track_d, year, lats, lons, write, out_path):
    df_1 = pd.DataFrame(track_d).T
    ids = [k for k in track_d.keys()]
    if len(ids) == 0:
        df_out = pd.DataFrame({'ID': [], 'cx_p': [], 'cy_p': [], 'c_p (Pa)': [], 
                               'cx_c': [], 'cy_c': []})
    else:
        df_out = pd.DataFrame({
            'ID': ids,
            'cx_p': df_1['cx_p'].to_numpy(),
            'cy_p': df_1['cy_p'].to_numpy(),
            'c_p (Pa)': df_1['c_p (Pa)'].to_numpy(),
            'cx_c': df_1['cx_c'].to_numpy(),
            'cy_c': df_1['cy_c'].to_numpy()
        })

    min_pressure = []
    total_distance_p = []
    total_distance_c = []
    lifetime = []
    time_min_p = []

    ny, nx = lats.shape  # Get grid shape once

    for i in range(len(ids)):
        if df_out['cx_p'][i][-1] == 'f':  # completed track
            lifetime_i = dt * (len(df_out['cx_p'][i]) - 2)
            time_min_p.append(dt * np.argmin(df_out['c_p (Pa)'][i]))
            min_pressure.append(np.min(df_out['c_p (Pa)'][i][:-1]))
            xpci = df_out['cx_p'][i][:-1]
            ypci = df_out['cy_p'][i][:-1]
            xcci = df_out['cx_c'][i][:-1]
            ycci = df_out['cy_c'][i][:-1]
        else:  # incomplete track
            lifetime_i = dt * (len(df_out['cx_p'][i]) - 1)
            time_min_p.append(dt * np.argmin(df_out['c_p (Pa)'][i]))
            min_pressure.append(np.min(df_out['c_p (Pa)'][i]))
            xpci = df_out['cx_p'][i]
            ypci = df_out['cy_p'][i]
            xcci = df_out['cx_c'][i]
            ycci = df_out['cy_c'][i]

        dpi = 0
        dci = 0
        for j in range(len(xpci) - 1):
            # --- Pressure Center (p) ---
            x1, y1 = xpci[j], ypci[j]
            x2, y2 = xpci[j+1], ypci[j+1]

            if not (0 <= x1 < nx and 0 <= y1 < ny and 0 <= x2 < nx and 0 <= y2 < ny):
                print(f"[P] Skipping out-of-bounds at j={j}: ({x1},{y1}) → ({x2},{y2})")
                continue

            lonp1, latp1 = lat_lon(x1, y1, lons, lats)
            lonp2, latp2 = lat_lon(x2, y2, lons, lats)

            if np.any(np.isnan([lonp1, latp1, lonp2, latp2])):
                print(f"[P] Skipping NaNs at j={j}")
                continue

            dp = phys_dist(latp2, lonp2, latp1, lonp1)
            dpi += dp
            print(f"[P] Step {j}: Distance = {float(dp):.2f} km")

            # --- Spatial Center (c) ---
            x1c, y1c = xcci[j], ycci[j]
            x2c, y2c = xcci[j+1], ycci[j+1]

            if not (0 <= x1c < nx and 0 <= y1c < ny and 0 <= x2c < nx and 0 <= y2c < ny):
                print(f"[C] Skipping out-of-bounds at j={j}: ({x1c},{y1c}) → ({x2c},{y2c})")
                continue

            lonc1, latc1 = lat_lon(x1c, y1c, lons, lats)
            lonc2, latc2 = lat_lon(x2c, y2c, lons, lats)

            if np.any(np.isnan([lonc1, latc1, lonc2, latc2])):
                print(f"[C] Skipping NaNs at j={j}")
                continue

            dc = phys_dist(latc2, lonc2, latc1, lonc1)
            dci += dc
            print(f"[C] Step {j}: Distance = {float(dc):.2f} km")

        lifetime.append(lifetime_i)
        total_distance_p.append(dpi)
        total_distance_c.append(dci)

    df_out['min_p (Pa)'] = min_pressure
    df_out['tot_dist_p (km)'] = total_distance_p
    df_out['tot_dist_c (km)'] = total_distance_c
    df_out['lifetime (h)'] = lifetime
    df_out['speed_p (km/h)'] = np.array(total_distance_p) / np.array(lifetime)
    df_out['speed_c (km/h)'] = np.array(total_distance_c) / np.array(lifetime)
    df_out['time_min_p (h)'] = time_min_p

    # Filter for storms passing through USEC region
    df_out = filter_USEC(df_out, lats, lons)
    print("=== USEC Filter Results ===")
    print(df_out[['ID', 'passed_USEC']])
    print("Num storms passing USEC:", (df_out['passed_USEC'] > 0).sum())

    # Optional write before distance filtering (to compare raw and filtered)
    if write and len(ids) != 0:
        out_path = out_path + str(year) + '_stormtracks_unfiltered.csv'
        df_out.to_csv(out_path)

    return df_out


# takes storm tracking dataframe and drops storms that don't pass USEC box
def USEC_only(df_out, year, out_path):
    df_out['passed_USEC'] = df_out['passed_USEC'].astype(int)
    df_out = df_out.drop(df_out[df_out['passed_USEC'] == 0].index)
    df_out = df_out.reset_index(drop=True)
    out_path = out_path + str(year) + '_stormtracks_filtered.csv'
    df_out.to_csv(out_path)
    return df_out

# =============================================================================

# Main storm tracking functions

# returns a list of gridcoords and their pressure values that have
# pressure <= filter_mslp and are the minimum pressure within detect_box_rad
def p_filter_box(data_arr, lons, lats): #, lati, loni):
    where_pressure = np.argwhere(data_arr <= filter_mslp)
    pvals = []
    gridcoords = []

    for gridcoord in where_pressure:
        gridy = gridcoord[0]
        gridx = gridcoord[1]
        pval_c = data_arr[gridy, gridx]
#        lon, lat = lat_lon(gridx, gridy, lons, lats)
#        x_rad = grid_rad_lon(gridy, detect_box_rad)
        _, lat_val = lat_lon(gridx, gridy, lons, lats)
        x_rad = grid_rad_lon(lat_val, detect_box_rad)
        xrange, yrange, borders = grow_box(gridx, gridy, x_rad, grid_rad_lat, data_arr)
        xrange = np.int64(xrange)
        xrange = xrange.flatten()
        yrange = yrange.flatten()
        # Create meshgrid for box indices to extract the 2D region from data_arr
        X, Y = np.meshgrid(xrange, yrange)
        # Use the meshgrid to extract the data from the 2D array
        box = data_arr[Y, X]  # Correct indexing of the data array
        # Find the index of the minimum pressure value in the box
        min_index = np.unravel_index(np.argmin(box), box.shape)

        # Compare pval_c with the minimum value in the box
        if (min_index == (gridy - yrange[0], gridx - xrange[0])):
            pvals.append(pval_c)
            gridcoords.append(gridcoord)

    # checkpoint: print the number of candidate pvals
    print('Number of candidate storm centers: ' + str(len(pvals)))

    return np.array(gridcoords), np.array(pvals)


# Search for and identify storm candidates via MSLP objects
def get_objects_mslp(t, data, lons, lats, times, restrictLat):
    time_1 = time.strftime("%m/%d/%Y, %H:%M:%S", time.localtime())
    t1 = time.time()
    print("Start: ", time_1)

    # Search for candidate grid points that have pressure < 1010 hPa and are 
    # the lowest within 1000 km
    data_arr = data[t]
    gridcoords, pvals = p_filter_box(data_arr, lons, lats)

    if len(pvals) == 0:
        df_pfilt = pd.DataFrame({'cx_p': [], 'cy_p': [], 'c_p (Pa)': []})
    else:
        df_pfilt = pd.DataFrame({'cx_p': gridcoords[:, 1], 'cy_p': gridcoords[:, 0], 
                                 'c_p (Pa)': pvals})
        df_pfilt = df_pfilt.sort_values(by=['c_p (Pa)']).reset_index(drop=True)

    # Checkpoint: print df_pfilt
    print(df_pfilt)

    # Construct feature objects consisting of all grid points within p_contour 
    # of candidate grid point
    # Also record spatial centers and gives an ID to each candidate
    xccs = []
    yccs = []
    ids = []

    for i in range(len(df_pfilt)):
        xi = int(df_pfilt['cx_p'][i])
        yi = int(df_pfilt['cy_p'][i])
        coords = np.nonzero(ski.segmentation.flood(data_arr, (yi, xi), 
                                                   tolerance=p_contour))

        # Spatial center is an average of the coordinates of the 200 Pa closed 
        # contour, weighted by 1/pressure
        w_sum = 0
        w_x = []
        w_y = []
        for j in range(len(coords[0])):
            pij = data_arr[coords[0][j], coords[1][j]]
            wij = 1 / pij
            w_sum += wij
            w_y.append(wij * coords[0][j])
            w_x.append(wij * coords[1][j])
        xccs.append(round(np.sum(w_x) / w_sum))
        yccs.append(round(np.sum(w_y) / w_sum))
        
        # Check if times[t] is a numpy.datetime64 object
        if isinstance(times[t], np.datetime64):
            # Convert numpy.datetime64 to a Python datetime object
            python_datetime = pd.to_datetime(times[t].astype(str))
            # Now apply strftime to format it as a string
            time_str = python_datetime.strftime("%m%d%Y-%H")
        else:
            # If it's already a string, use it directly
            time_str = str(times[t])

        # Append the ID using the string format
        ids.append(time_str + '-' + str(xi) + '-' + str(yi))

    df_pfilt['cx_c'] = xccs
    df_pfilt['cy_c'] = yccs
    df_pfilt['ID'] = ids

    # Filter if pressure center and spatial center are too far away from each other
    xpcs = df_pfilt['cx_p']
    ypcs = df_pfilt['cy_p']
    xccs = df_pfilt['cx_c']
    yccs = df_pfilt['cy_c']
    pc_diff = []
    for i in range(len(xccs)):
        lonc, latc = lat_lon(xccs[i], yccs[i], lons, lats)
        print('lonc: ' + str(lonc) + '     latc: ' + str(latc))
        lonp, latp = lat_lon(xpcs[i], ypcs[i], lons, lats)
        print('lonp: ' + str(lonp) + '     latp: ' + str(latp))
        pc_diff.append(phys_dist(latc, lonc, latp, lonp))
    df_pfilt['pc_diff'] = pc_diff
    
    # Checkpoint: print pc_diff
    print('pc_diff: ' + str(pc_diff))
    
    df_pfilt = df_pfilt.drop(df_pfilt[df_pfilt['pc_diff'] >= pc_center_diff].index)
    df_pfilt = df_pfilt.reset_index(drop=True)

    # Checkpoint: print df_pfilt
    print(df_pfilt)

    time_2 = time.strftime("%m/%d/%Y, %H:%M:%S", time.localtime())
    t2 = time.time()
    print("End: ", time_2)
    print(f"Time elapsed: {t2 - t1} s")
    
    return df_pfilt


# input: loaded data for 1 year (sliced beforehand); applies full tracking
def full_run_1yr_data_loaded(year, datai, timesi, latsi, lonsi):
    # array of timesteps to run
    tarri = np.arange(0, len(timesi)) 
    print('time array: ', tarri*dt) 

    # initialize for start time
    df_t0 = get_objects_mslp(tarri[0], datai, lonsi, latsi, timesi, False)
    storm_tracks = {}
    dlist0 = key_val(df_t0)
    for i in range(len(df_t0)):
        storm_tracks[df_t0['ID'][i]] = dlist0[i]
    print(storm_tracks)
    print("\n")

    # timestepping
    for t in tarri[1:]:
        print(f't_hour: {t}')
        df_ti = get_objects_mslp(t, datai, lonsi, latsi, timesi, False)
        storms_not_ended = []
        for i in range(len(df_ti)):
            xip = df_ti['cx_p'][i]
            yip = df_ti['cy_p'][i]
            pi = df_ti['c_p (Pa)'][i]
            xi = df_ti['cx_c'][i]
            yi = df_ti['cy_c'][i]
            idi = df_ti['ID'][i]
            added = 0
            for k,v in storm_tracks.items():
                xvp = v['cx_p'][-1]
                if xvp == 'f':
                    continue
                yvp = v['cy_p'][-1]
                pv = v['c_p (Pa)'][-1]
                xv = v['cx_c'][-1]
                yv = v['cy_c'][-1]
                lonic, latic = lat_lon(xi, yi, lonsi, latsi)
                lonv, latv = lat_lon(xv, yv, lonsi, latsi)
                dist = phys_dist(latic, lonic, latv, lonv)
                if (dist <= search_rad) and (np.abs(pv-pi) <= slp_rad):
                    v['cx_p'].append(xip)
                    v['cy_p'].append(yip)
                    v['c_p (Pa)'].append(pi)
                    v['cx_c'].append(xi)
                    v['cy_c'].append(yi)
                    added = 1
                    storms_not_ended.append(k)
                    break

            # if candidate is not added to an existing track, then add new 
            # entry to storm_tracks dictionary
            if (added == 0): 
                val = {}
                val['cx_p'] = [xip]
                val['cy_p'] = [yip]
                val['c_p (Pa)'] = [pi]
                val['cx_c'] = [xi]
                val['cy_c'] = [yi]
                storm_tracks[idi] = val
                storms_not_ended.append(idi)

        # if track is complete, append 'f'; or delete the track if needed
        storm_to_del = []
        for k,v in storm_tracks.items():
            if (k not in storms_not_ended) and (v['cx_p'][-1] != 'f'):
                v['cx_p'].append('f')
                v['cy_p'].append('f')
                v['c_p (Pa)'].append('f')
                v['cx_c'].append('f')
                v['cy_c'].append('f')
                if (len(v['cx_p']) < t_req):
                    storm_to_del.append(k)
        for key in storm_to_del:
            del storm_tracks[key]
        # Yield the updated storm tracks at each timestep
        yield storm_tracks  # Returns the current state of storm tracks
        print("\n")

    return storm_tracks

# =============================================================================

# Global variables

# initial tracking parameters
detect_box_rad = 1000 # km
grid_rad_lat = round(detect_box_rad/(111.3 * 0.11))  # CRCM5 has 0.11 deg (~12 km) spatial resolution
filter_mslp = 101000 # in Pa
p_contour = 200 # Pa, closed low pressure feature

# spatial limits for nor'easter formation (if using filter_gen)
max_longitude = -95 # deg lon
max_latitude = 20 # deg lat

# limit for post-processing psl filter
minSLP = 98000 # in Pa

# lat/lon bounding box for US east coast
usne_N = 45 # deg lat
usne_S = 33 # deg lat
usne_E = -70 # deg lon
usne_W = -80 # deg lon

# north pole location in CRCM5 rotated coordinate system
north_pole = [83.0, 42.5]
pole_lon = 83.0
pole_lat = 42.5

# timestepping parameters
dt = 6 # h
filter_lifetime = 24 # in h
t_req = int(filter_lifetime/dt + 2)
search_rad = 125*dt # km, assumes NHemi ETCs can reach speeds of 125 km/h, multiply by length of timestep
slp_rad = 10000 # Pa
filter_distance = 1000 # in km
pc_center_diff = 500 # km

# =============================================================================

'''

# Usage example

# Now process the tracking for a single storm season (Sept [Year] to April [Year+1]), 
# keyed by end year
tracks_data = full_run_1yr_data_loaded({year+1}, storm_season.psl.values, 
                                       storm_season.time.values, 
                                       storm_season.lat.values, 
                                       storm_season.lon.values)

# Convert tracks_data generator into a dictionary for use in post-processing
final_tracks = None
for tracks in tracks_data:
    final_tracks = tracks  # Keep the latest iteration (final version of the storm tracks)

# Post-processing pipeline
tracks_df = storm_track_out(final_tracks, {year+1}, storm_season.lat.values, 
                            storm_season.lon.values, 1, unfiltered_out_path)

# Apply the post-processing filters
filter_tracks_df = SLP_filter(tracks_df, minSLP)
filter_tracks_df = filter_dist(filter_tracks_df)
filter_tracks_df = filter_time(filter_tracks_df)
filter_tracks_df = USEC_only(filter_tracks_df, {year+1}, filtered_out_path)

'''

# =============================================================================



