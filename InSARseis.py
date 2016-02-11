
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

    helptext = 'Plot only precalculated NumPy file'
    parser.add_argument('-p', '--plot_only', help=helptext, default=False,
                        action='store_true')

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

def plot_result(src_file):
    import instaseis
    import matplotlib.pyplot as plt

    npy_file = '%s.npy' % (os.path.splitext(src_file)[0])
    all_data = np.load(npy_file)

    source = instaseis.FiniteSource.from_usgs_param_file(src_file)
    pntsrc_lats = np.zeros(source.npointsources)
    pntsrc_lons = np.zeros(source.npointsources)
    pntsrc_depth = np.zeros(source.npointsources)
    pntsrc_M0 = np.zeros(source.npointsources)
    for i in range(0, source.npointsources):
        pntsrc_depth[i] = source.pointsources[i].depth_in_m
        pntsrc_lats[i] = source.pointsources[i].latitude
        pntsrc_lons[i] = source.pointsources[i].longitude
        pntsrc_M0[i] = source.pointsources[i].M0

    lats = all_data[0]
    lons = all_data[1]
    offset_E = all_data[2]
    offset_N = all_data[3]
    offset_Z = all_data[4]
    
    vmin = -np.around(np.max(np.abs((offset_E, offset_N, offset_Z)))/5, 1)
    vmax = np.round(np.max(np.abs((offset_E, offset_N, offset_Z)))/5, 1)

    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(221)
    #h1 = ax1.pcolormesh(lons, lats, offset_E, vmin=vmin, vmax=vmax, cmap='seismic')
    h1 = plot_on_map(lons, lats, offset_E, vmin, vmax, ax=ax1)
    ax1.set_title('E/W displacement')

    ax2 = fig.add_subplot(222)
    #h2 = ax2.pcolormesh(lons, lats, offset_N, vmin=vmin, vmax=vmax, cmap='seismic')
    h2 = plot_on_map(lons, lats, offset_N, vmin, vmax, ax=ax2)
    ax2.set_title('N/S displacement')

    ax3 = fig.add_subplot(223)
    #h3 = ax3.pcolormesh(lons, lats, offset_Z, vmin=vmin, vmax=vmax, cmap='seismic')
    h3 = plot_on_map(lons, lats, offset_Z, vmin, vmax, ax=ax3)
    ax3.set_title('vertical displacement')

    ax4 = fig.add_subplot(224)
    #h4 = ax4.scatter(y=pntsrc_lats, x=pntsrc_lons, c=pntsrc_depth/1e3, linewidths=0, cmap='GnBu')
    h4 = plot_fault_on_map(pntsrc_lons, pntsrc_lats, lons, lats, 
                           pntsrc_M0/pntsrc_M0.max()*40., pntsrc_depth/1e3, 
                           vmin, vmax, ax=ax4)
    ax4.set_xlim(ax3.get_xlim())
    ax4.set_ylim(ax3.get_ylim())
    ax4.set_title('Fault')

    for ax in fig.axes:
        ax.set_ylabel('Latitude')
        ax.set_xlabel('Longitude')
    ax_cbar_disp = fig.add_axes([0.87, 0.536, 0.035, 0.365])
    cbar_disp = fig.colorbar(h1, cax=ax_cbar_disp)
    cbar_disp.set_label('displacement / m')
    
    ax_cbar_depth = fig.add_axes([0.87, 0.1, 0.035, 0.365])
    cbar_depth = fig.colorbar(h4, cax=ax_cbar_depth)
    cbar_depth.set_label('depth / km')

    fig.subplots_adjust(right=0.83)

    
    png_file = '%s.png' % (os.path.splitext(npy_file)[0])
    print 'Saving to file %s' % png_file
    
    fig.savefig(png_file, frameon=False)
    plt.show()
    
def plot_on_map(lons, lats, data, vmin, vmax, ax):
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='stere', 
                lat_0=np.mean(lats), lon_0=np.mean(lons), 
                llcrnrlat=lats.min(), urcrnrlat=lats.max(),
                llcrnrlon=lons.min(), urcrnrlon=lons.max(),
                resolution='i', ax=ax, area_thresh=1000.)
    m.drawcoastlines(linewidth=0.35)
    m.drawcountries(linewidth=0.25)
    parallels = np.arange(-90.,90,5.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = np.arange(0.,360.,5.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    x, y = m(lons, lats)

    h = m.pcolormesh(x, y, data=data, vmin=vmin, vmax=vmax, 
                     cmap='seismic', ax=ax)

    return h

def plot_fault_on_map(lons_fault, lats_fault, lons_mesh, lats_mesh,
                      data_1, data_2, vmin, vmax, ax):
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='stere', 
                lat_0=np.mean(lats_mesh), lon_0=np.mean(lons_mesh), 
                llcrnrlat=lats_mesh.min(), urcrnrlat=lats_mesh.max(),
                llcrnrlon=lons_mesh.min(), urcrnrlon=lons_mesh.max(),
                resolution='i', ax=ax, area_thresh=1000.)
    m.drawcoastlines(linewidth=0.35)
    m.drawcountries(linewidth=0.25)

    # draw parallels
    parallels = np.arange(-90., 90., 5.)
    m.drawparallels(parallels, labels=[1,0,0,0], fontsize=10)

    # draw meridians
    meridians = np.arange(0., 360., 5.)
    m.drawmeridians(meridians, labels=[0,0,0,1], fontsize=10)

    x, y = m(lons_fault, lats_fault)
    h = m.scatter(x, y, s=data_1, c=data_2, linewidths=0, cmap='GnBu')
        
    return h
    

if __name__ == '__main__':
    args = define_argument_parser()

    if not args.plot_only:
      all_data = calc_static_offset(args.source_file_name, args.database_path)

    plot_result(args.source_file_name)


