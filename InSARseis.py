
# coding: utf-8

# In[74]:


import matplotlib.pyplot as plt
import numpy as np
import os
from progressbar import ProgressBar, ETA, Bar, Percentage
import argparse

def define_argument_parser():
    helptext = 'Calculate Static offset around an earthquake.'
    formatter_class = argparse.RawTextHelpFormatter
    parser = argparse.ArgumentParser(description=helptext,
                                     formatter_class=formatter_class)

    helptext = 'Source file in USGS Param format'
    parser.add_argument('-s', '--source_file_name', help=helptext, required=True)

    helptext = 'Path to Instaseis database'
    parser.add_argument('-d', '--database_path', help=helptext, required=True)

    args = parser.parse_args()
    return args

def calc_static_offset(src_file, db_dir, nx=400, ny=400):
    import obspy
    import instaseis

    db = instaseis.open_db(db_dir)
    print db

    source = instaseis.FiniteSource.from_usgs_param_file(src_file, dt=db.info.dt)
    print source

    pntsrc_lats = np.zeros(source.npointsources)
    pntsrc_lons = np.zeros(source.npointsources)
    pntsrc_depth = np.zeros(source.npointsources)
    
    for i in range(0, source.npointsources):
        pntsrc_depth[i] = source.pointsources[i].depth_in_m
        pntsrc_lats[i] = source.pointsources[i].latitude
        pntsrc_lons[i] = source.pointsources[i].longitude

    x_extent = source.max_longitude - source.min_longitude
    y_extent = source.max_latitude - source.min_latitude
    extent = np.max((x_extent, y_extent))

    latmin = source.hypocenter_latitude - extent * 1.5
    lonmin = source.hypocenter_longitude - extent * 1.5
    latmax = source.hypocenter_latitude + extent * 1.5
    lonmax = source.hypocenter_longitude + extent * 1.5

    print 'Window: ', latmin, latmax, lonmin, lonmax

    lats = []
    lons = []
    offset = []
    lats, lons = np.meshgrid(np.linspace(latmin, latmax, nx), 
                             np.linspace(lonmin, lonmax, ny))
    offset_E = np.zeros_like(lats)
    offset_N = np.zeros_like(lats)
    offset_Z = np.zeros_like(lats)

    # Define progress bar
    widgets = ['Calculating displacement: ', Percentage(), ' ',
               Bar(), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=nx).start()

    for ix in range(0, lats.shape[0]):
        for iy in range(0, lats.shape[1]):
            npts = db.info.npts
            rec = instaseis.Receiver(latitude=lats[ix, iy], longitude=lons[ix, iy])
            st = db.get_seismograms_finite_source(sources=source, receiver=rec, components='ENZ')
            offset_E[ix, iy] = np.mean(st.select(channel='LXE')[0].data[npts*0.5:npts])
            offset_N[ix, iy] = np.mean(st.select(channel='LXN')[0].data[npts*0.5:npts])
            offset_Z[ix, iy] = np.mean(st.select(channel='LXZ')[0].data[npts*0.5:npts])        

        # Update Progress Bar
        pbar.update(ix)

        all_data = [lats, lons, offset_E, offset_N, offset_Z]
        fnam = '%s.npy' % (os.path.splitext(src_file)[0])
        np.save(fnam, all_data, allow_pickle=True)

    return all_data

def plot_result(npy_file):

    all_data = np.load(npy_file)
    lats = all_data[0]
    lons = all_data[1]
    offset_E = all_data[2]
    offset_N = all_data[3]
    offset_Z = all_data[4]
    
    vmin = -np.around(np.max(np.abs((offset_E, offset_N, offset_Z)))/10, 1)
    vmax = np.round(np.max(np.abs((offset_E, offset_N, offset_Z)))/10, 1)

    fig = plt.figure(figsize=(15,15))
    ax1 = plt.subplot(221)
    h1 = ax1.pcolormesh(lons, lats, offset_E, vmin=vmin, vmax=vmax, cmap='seismic')
    ax1.set_title('E/W displacement')

    ax2 = plt.subplot(222)
    h2 = ax2.pcolormesh(lons, lats, offset_N, vmin=vmin, vmax=vmax, cmap='seismic')
    ax2.set_title('N/S displacement')

    ax3 = plt.subplot(223)
    h3 = ax3.pcolormesh(lons, lats, offset_Z, vmin=vmin, vmax=vmax, cmap='seismic')
    ax3.set_title('vertical displacement')

    ax4 = plt.subplot(224)
    h4 = ax4.scatter(y=pntsrc_lats, x=pntsrc_lons, c=pntsrc_depth, linewidths=0)
    ax4.set_xlim(ax3.get_xlim())
    ax4.set_ylim(ax3.get_ylim())
    ax4.set_title('Fault')

    for ax in fig.axes:
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Longitude')
    ax_cbar = fig.add_axes([0.87, 0.15, 0.035, 0.7])
    cbar = fig.colorbar(h1, cax=ax_cbar)
    fig.subplots_adjust(right=0.83)
    cbar.set_label('displacement / m')
    
    fnam = '%s.png' % (os.path.splitext(npy_file)[0])
    
    fig.savefig(fnam)

if __name__ == '__main__':
    args = define_argument_parser()

    all_data = calc_static_offset(args.source_file_name, args.database_path)


