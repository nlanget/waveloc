import os
import h5py
import logging
import numpy as np
import matplotlib.pyplot as plt
from waveloc.filters import smooth
from waveloc.NllGridLib import read_hdr_file, read_stations_file
from waveloc.options import WavelocOptions
from synth_migration import generateSyntheticDirac
from waveloc.locations_trigger import trigger_locations_inner

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s : %(asctime)s : %(message)s')


def setUp():

    logging.info('Setting up synthetic test case generation...')

    # set up default parameters
    wo = WavelocOptions()

    # set base path to $WAVELOC_PATH
    wo.verify_base_path()

    # set options for synthetic test
    wo.opdict['time'] = True
    wo.opdict['verbose'] = False

    wo.opdict['outdir'] = 'EXAMPLE_Dirac'
    wo.opdict['time_grid'] = 'Slow_len.100m.P'
    wo.opdict['load_ttimes_buf'] = True
    wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    wo.opdict['stations'] = 'coord_stations_piton'

    wo.opdict['syn_filename'] = 'test_grid4D_hires.hdf5'
    wo.opdict['syn_amplitude'] = [1.0]
    wo.opdict['syn_datalength'] = 20.0
    wo.opdict['syn_samplefreq'] = 100.0
    wo.opdict['syn_kwidth'] = 0.1
    wo.opdict['syn_otime'] = [6.0]

    # place default point at center of grid
    base_path = wo.opdict['base_path']
    search_grid_filename = os.path.join(base_path, 'lib',
                                        wo.opdict['search_grid'])
    grid_info = read_hdr_file(search_grid_filename)
    wo.opdict['syn_ix'] = [grid_info['nx']/2]
    wo.opdict['syn_iy'] = [grid_info['ny']/2]
    wo.opdict['syn_iz'] = [grid_info['nz']/2]

    # sanity check
    wo.verify_synthetic_options()

    return (wo, grid_info)


def doPointTest(wo, loclevel=10.0):

    logging.info('Doing synthetic test for point (%d,%d,%d)...' %
                 (wo.opdict['syn_ix'][0], wo.opdict['syn_iy'][0],
                  wo.opdict['syn_iz'][0]))

    # do the migration
    test_info = generateSyntheticDirac(wo.opdict)
    logging.info(test_info)

    # retrieve output info
    stack_filename = test_info['stack_file']
    nx, ny, nz, nt = test_info['grid_shape']
    dx, dy, dz, dt = test_info['grid_spacing']
    x_orig, y_orig, z_orig = test_info['grid_orig']
    ix_true, iy_true, iz_true, it_true = test_info['true_indexes']
    stack_start_time = test_info['start_time']

    # extract the max stacks
    f_stack = h5py.File(stack_filename, 'r')
    max_val = f_stack['max_val']
    max_x = f_stack['max_x']
    max_y = f_stack['max_y']
    max_z = f_stack['max_z']

    # smooth the max_val
    max_val_smoothed = smooth(max_val[:])

    # launch a location trivver
    locs = trigger_locations_inner(max_val_smoothed, max_x, max_y, max_z,
                                   loclevel, loclevel, stack_start_time, dt)

    f_stack.close()
    return test_info, locs


def analyseLocs(locs, wo, test_info):

    dx, dy, dz, dt = test_info['grid_spacing']
    x_orig, y_orig, z_orig = test_info['grid_orig']
    ix_true, iy_true, iz_true, it_true = test_info['true_indexes']

    n_locs = len(locs)

    #fig = plt.figure()
    #fig.set_facecolor('white')
    #for l,key in enumerate(sorted(test_info['data'])):
    #  ax = fig.add_subplot(len(sorted(test_info['data'])),1,l+1)
    #  ax.plot(test_info['data'][key],'k')
    #  ax.text(0.2,0.5,key)
    #  ax.set_axis_off()
    #plt.show()


    loc_dt_list = [aloc['o_time'] - wo.opdict['syn_otime'] for aloc in locs]
    loc_dist_list = [np.sqrt((ix_true*dx+x_orig - aloc['x_mean'])**2 +
                             (iy_true*dy+y_orig - aloc['y_mean'])**2 +
                             (iz_true*dy+z_orig - aloc['z_mean'])**2)
                     for aloc in locs]

    # This is a dirac test, but we may have more than one loc
    # (secondary maxima) so for safety, pull out best loc
    if n_locs > 0:
        imax = np.argmax([loc['max_trig'] for loc in locs])
        trig_loc = locs[imax]
        loc_dt = loc_dt_list[imax]
        loc_dist = loc_dist_list[imax]
    else:
        trig_loc = None
        loc_dist = None
        loc_dt = None

    if n_locs > 1:
      for loc in locs:
        plotSyntheticResult(wo,test_info,loc)
      plt.show()


    return n_locs, loc_dist, loc_dt, trig_loc


def doResolutionTest(wo, grid_info, filename, loclevel=10.0,
                     decimation=(1, 1, 1)):

    dec_x, dec_y, dec_z = decimation

    # get information from search grid
    nx = grid_info['nx']
    ny = grid_info['ny']
    nz = grid_info['nz']

    # get decimated grid
    nx_dec = int(nx/dec_x)
    ny_dec = int(ny/dec_y)
    nz_dec = int(nz/dec_z)
    nb_dec = nx_dec*ny_dec*nz_dec

    # create info for decimated grid
    dec_grid_info = {}
    dec_grid_info['nx'] = nx_dec
    dec_grid_info['ny'] = ny_dec
    dec_grid_info['nz'] = nz_dec
    dec_grid_info['dx'] = grid_info['dx']*dec_x
    dec_grid_info['dy'] = grid_info['dy']*dec_y
    dec_grid_info['dz'] = grid_info['dz']*dec_z
    dec_grid_info['x_orig'] = grid_info['x_orig']
    dec_grid_info['y_orig'] = grid_info['y_orig']
    dec_grid_info['z_orig'] = grid_info['z_orig']

    # set up arrays to store values
    dist_grid = np.empty(nb_dec, dtype='float32')
    nloc_grid = np.empty(nb_dec, dtype='int16')
    dt_grid = np.empty(nb_dec, dtype='float32')
    dist_grid[...] = -999.0
    nloc_grid[...] = -999
    dt_grid[...] = -999.0

    # iterate over points
    for ib_dec in xrange(nb_dec):
        ix_dec, iy_dec, iz_dec = np.unravel_index(ib_dec,
                                                  (nx_dec, ny_dec, nz_dec))
        # reconstruct the indexes in the original grid
        ix = ix_dec*dec_x
        iy = iy_dec*dec_y
        iz = iz_dec*dec_z

        wo.opdict['syn_ix'] = [ix]
        wo.opdict['syn_iy'] = [iy]
        wo.opdict['syn_iz'] = [iz]

        # do synthetic test for this point
        test_info, locs = doPointTest(wo, loclevel=10.0)

        # do analysis for this point
        n_locs, best_dist, best_dt, trig_loc = analyseLocs(locs, wo, test_info)
        print ix, iy, iz, n_locs, best_dist, best_dt

        # save into array
        dist_grid[ib_dec] = best_dist
        nloc_grid[ib_dec] = n_locs
        dt_grid[ib_dec] = best_dt

    # write HDF5 file
    resgrid_filename = os.path.join(wo.opdict['base_path'], 'out',
                                    wo.opdict['outdir'], 'grid', filename)
    f = h5py.File(resgrid_filename, 'w')
    f_dist = f.create_dataset('dist_grid', data=dist_grid)
    f_nloc = f.create_dataset('nloc_grid', data=nloc_grid)
    f_dt = f.create_dataset('dt_grid', data=dt_grid)
    for key, value in dec_grid_info.iteritems():
        f_dist.attrs[key] = value
        f_nloc.attrs[key] = value
        f_dt.attrs[key] = value
    f.close()

    return resgrid_filename


def plotSyntheticResult(wo,test_info,loc):

            # The first figure shows the stacking values for the true origin time.
            # The second figure shows the stacking values for the recovered origin time.
            # The third plot shows the kurtosis traces alignment.
            # yellow star = true location ; red star = Waveloc's location

            import matplotlib.pyplot as plt

            nx, ny, nz, nt = test_info['grid_shape']
            dx, dy, dz, dt = test_info['grid_spacing']
            x_orig, y_orig, z_orig = test_info['grid_orig']
            ix_true, iy_true, iz_true, it_true = test_info['true_indexes']

            stack_filename = test_info['stack_file']
            f_stack = h5py.File(stack_filename,'r')
            if 'max_val_smooth' in f_stack:
              max_val = f_stack['max_val_smooth']
            else:
              max_val = f_stack['max_val']
            max_val = f_stack['max_val']

            ixx = wo.opdict['syn_ix'][0]
            iyy = wo.opdict['syn_iy'][0]
            izz = wo.opdict['syn_iz'][0]

            # Waveloc location
            it_found = wo.opdict['syn_samplefreq'] * loc['o_time']
            ix_found = int(round((loc['x_mean']-x_orig)*1./dx))
            iy_found = int(round((loc['y_mean']-y_orig)*1./dy))
            iz_found = int(round((loc['z_mean']-z_orig)*1./dz))

            from plot_mpl import plotDiracTest
            print "True location", test_info['true_indexes']
            fig_name = os.path.join(wo.opdict['base_path'], 'out',wo.opdict['outdir'],'fig')
            wo.opdict['syn_filename'] = 'loc1_x%d_y%d_z%d'%(ixx,iyy,izz)
            plotDiracTest(test_info,fig_name,2,ixx,iyy,izz,ix_found,iy_found,iz_found)

            test_info['true_indexes'] = (ix_found, iy_found, iz_found, it_found)
            print "Waveloc location",test_info['true_indexes']
            fig_name = os.path.join(wo.opdict['base_path'], 'out',wo.opdict['outdir'],'fig')
            wo.opdict['syn_filename'] = 'loc2_x%d_y%d_z%d'%(ixx,iyy,izz)
            plotDiracTest(test_info,fig_name,2,ixx,iyy,izz,ix_found,iy_found,iz_found)

            # True location
            wix = np.argmin(np.abs(np.array(wo.opdict['syn_ix'])*dx+x_orig - np.arange(x_orig,x_orig+nx*dx,dx)))
            wiy = np.argmin(np.abs(np.array(wo.opdict['syn_iy'])*dy+y_orig - np.arange(y_orig,y_orig+ny*dy,dy)))
            wiz = np.argmin(np.abs(np.array(wo.opdict['syn_iz'])*dz+z_orig - np.arange(z_orig,z_orig+nz*dz,dz)))
            wib = wix*ny*nz + wiy*nz + wiz

            from hdf5_grids import get_interpolated_time_grids
            time_grids = get_interpolated_time_grids(wo.opdict)
            ttimes_true = {}
            for sta in sorted(test_info['data']):
              if time_grids.has_key(sta):
                ttimes_true[sta] = time_grids[sta].grid_data[wib]
              else:
                logging.info('Missing travel-time information for station %s. Ignoring station...'%sta)

            # Waveloc location
            wib = ix_found*ny*nz + iy_found*nz + iz_found
            ttimes_wl = {}
            for sta in sorted(test_info['data']):
              if time_grids.has_key(sta):
                ttimes_wl[sta] = time_grids[sta].grid_data[wib]
              else:
                logging.info('Missing travel-time information for station %s. Ignoring station...'%sta)


            fig = plt.figure(figsize=(12,8))
            fig.set_facecolor('white')
            n_traces = len(ttimes_true) + 2
            list_ind = [] 
            for iii,key in enumerate(sorted(ttimes_true)):
              delta_t = test_info['grid_spacing'][-1]
              t = np.arange(0,len(test_info['data'][key])*delta_t,delta_t)
              ax = fig.add_subplot(n_traces,3,3*iii+1)
              ax.set_axis_off()
              ax.plot(t,test_info['data'][key],'k')
              ax.text(0.2,0.5,key)

              ind = np.argmin(np.abs(t-ttimes_true[key]))
              t = np.arange(0,len(test_info['data'][key][ind:ind+len(max_val[:])])*delta_t,delta_t)
              ax = fig.add_subplot(n_traces,3,3*iii+2)
              ax.set_axis_off()
              if len(t) == len(test_info['data'][key][ind:ind+len(max_val[:])]):
                ax.plot(t,test_info['data'][key][ind:ind+len(max_val[:])],'k')
              else:
                ax.plot(t[:-1],test_info['data'][key][ind:ind+len(max_val[:])],'k')

              ind = np.argmin(np.abs(t-ttimes_wl[key]))
              t = np.arange(0,len(test_info['data'][key][ind:ind+len(max_val[:])])*delta_t,delta_t)
              ax = fig.add_subplot(n_traces,3,3*iii+3)
              ax.set_axis_off()
              if len(t) == len(test_info['data'][key][ind:ind+len(max_val[:])]):
                ax.plot(t,test_info['data'][key][ind:ind+len(max_val[:])],'k')
              else:
                ax.plot(t[:-1],test_info['data'][key][ind:ind+len(max_val[:])],'k')

            t = np.arange(0,len(max_val[:])*delta_t,delta_t)
            ind = np.argmin(np.abs(t-loc['o_time']))
            ax = fig.add_subplot(n_traces,3,3*(iii+1)+2)
            ax.set_axis_off()
            t = np.arange(0,len(max_val[ind:])*delta_t,delta_t)
            ax.plot(max_val[:],'k')

            ax = fig.add_subplot(n_traces,3,3*(iii+1)+3)
            ax.set_axis_off()
            ax.plot(max_val[:],'k')

            plt.figtext(0.45,0.93,'True location')
            plt.figtext(0.72,0.93,'Waveloc location')

            f_stack.close()


def plotResolutionTest(wo, hdf_filename, plot_filename):

    # read HDF5 file
    full_filename = os.path.join(wo.opdict['base_path'], 'out',
                                 wo.opdict['outdir'], 'grid', hdf_filename)
    f = h5py.File(full_filename, 'r')

    # grilles
    dist_grid_hdf = f['dist_grid']
    nloc_grid_hdf = f['nloc_grid']
    dt_grid_hdf = f['dt_grid']

    # attributs
    grid_info = dist_grid_hdf.attrs
    nx = grid_info['nx']
    ny = grid_info['ny']
    nz = grid_info['nz']
    dx = grid_info['dx']
    dy = grid_info['dy']
    dz = grid_info['dz']
    x_orig = grid_info['x_orig']
    y_orig = grid_info['y_orig']
    z_orig = grid_info['z_orig']

    # versions ram des grilles
    dist_grid = dist_grid_hdf[:].reshape(nx, ny, nz)
    nloc_grid = nloc_grid_hdf[:].reshape(nx, ny, nz)
    dt_grid = dt_grid_hdf[:].reshape(nx, ny, nz)

    # close HDF5 files
    f.close()

    # set up axes
    x = np.arange(nx)*dx
    y = np.arange(ny)*dy

    # read station file
    fname = os.path.join(wo.opdict['base_path'], 'lib', wo.opdict['stations'])
    stations = read_stations_file(fname)
    sta_x = np.empty(len(stations), dtype='float32')
    sta_y = np.empty(len(stations), dtype='float32')
    i = 0
    for key, value in stations.iteritems():
        if value['loc_type'] == 'XYZ':
            sta_x[i] = value['x']-x_orig
            sta_y[i] = value['y']-y_orig
        else:
            sta_x[i] = value['lon']-x_orig
            sta_y[i] = value['lat']-y_orig
        i = i+1

    # do plot
    vmin_dist = np.min(dist_grid)
    vmax_dist = np.max(dist_grid)
    vmin_nloc = np.min(nloc_grid)
    vmax_nloc = np.max(nloc_grid)
    vmin_dt = np.min(dt_grid)
    vmax_dt = np.max(dt_grid)

    # filename
    (root, ext) = os.path.splitext(plot_filename)

    col = plt.cm.jet
    # iterate over iz
    for iz in xrange(nz):
        plt.clf()
        plt.figure(figsize=(10, 4.5))

        # get depth
        zvalue = iz*dz+z_orig
        fname = root+'_%.2fkm' % zvalue+ext
        full_fname = os.path.join(wo.opdict['base_path'], 'out',
                                  wo.opdict['outdir'], 'fig', fname)

        xy_cut_dist = dist_grid[:, :, iz]
        xy_cut_nloc = nloc_grid[:, :, iz]
        xy_cut_dt = dt_grid[:, :, iz]
        xy_extent = [np.min(x), np.max(x), np.min(y), np.max(y)]

        p = plt.subplot(1, 3, 1)
        plt.imshow(xy_cut_dist.T, vmin=vmin_dist, vmax=vmax_dist,
                   origin='lower', interpolation='none',
                   extent=xy_extent, cmap=col)
        p.tick_params(labelsize=10)
        p.xaxis.set_ticks_position('bottom')
        plt.xlabel('x (km wrt ref)', size=10)
        plt.ylabel('y (km wrt ref)', size=10)
        plt.colorbar(orientation='horizontal', ticks=[0, 0.1, 0.2, 0.3, 0.4])
        plt.scatter(sta_x, sta_y)
        plt.figtext(0.08,0.71,'(a)')
        plt.title('Distance (km)')

        p = plt.subplot(1, 3, 2)
        plt.imshow(xy_cut_nloc.T, vmin=vmin_nloc, vmax=vmax_nloc,
                   origin='lower', interpolation='none',
                   extent=xy_extent, cmap=col)
        p.tick_params(labelsize=10)
        p.xaxis.set_ticks_position('bottom')
        plt.xlabel('x (km wrt ref)', size=10)
        plt.colorbar(orientation='horizontal', ticks=[0, 1, 2, 3])
        plt.scatter(sta_x, sta_y)
        plt.figtext(0.38,0.71,'(b)')
        plt.title('No of locs')

        p = plt.subplot(1, 3, 3)
        plt.imshow(xy_cut_dt.T, vmin=vmin_dt, vmax=vmax_dt,
                   origin='lower', interpolation='none',
                   extent=xy_extent, cmap=col)
        p.tick_params(labelsize=10)
        p.xaxis.set_ticks_position('bottom')
        plt.xlabel('x (km wrt ref)', size=10)
        plt.colorbar(orientation='horizontal', ticks=np.arange(0,vmax_dt,.025))
        plt.scatter(sta_x, sta_y)
        plt.figtext(0.65,0.71,'(c)')
        plt.title('Dt origin time (s)')

        plt.figtext(0.12,0.78,'z=%.2fkm'%zvalue)
        plt.savefig(full_fname)

if __name__ == '__main__':

    hdf_filename = 'waveloc_resolution.hdf5'
    plot_filename = 'waveloc_resolution.png'

    wo, grid_info = setUp()
    #doResolutionTest(wo, grid_info, hdf_filename, loclevel=10.0, decimation=(5, 5, 3))
    plotResolutionTest(wo, hdf_filename, plot_filename)
