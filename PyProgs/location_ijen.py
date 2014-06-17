import numpy as np
import pandas as pd
import os, sys, glob
import logging


def list_of_files():

  filename = '../../discrimination/lib/Ijen/Ijen_reclass_all_ter.csv'
  datapath = '/media/disk/CLASSIFICATION'
  list_sta = ['DAM','IJEN','KWUI','MLLR','POS','POSI','PSG','RAUN','TRWI']
  dest_path = '../data/Ijen'

  df = pd.read_csv(filename)
  df = df[df.Type=='VulkanikB']

  new_df = pd.DataFrame(index=df.Date.values,columns=list_sta)
  for sta in list_sta:
    station = sta.lower()
    for date in new_df.index:
      files = glob.glob(os.path.join(datapath,station,'*Z*%s_%s*'%(str(date)[:8],str(date)[8:])))
      if len(files) > 0:
        file = files[0]
        new_df[sta][new_df.index==date] = 1
        #os.symlink(file,os.path.join(dest_path,os.path.basename(file)))
      else:
        new_df[sta][new_df.index==date] = 0
    
  new_df.to_csv('../../discrimination/lib/Ijen/vb_sta.csv')


def ijen():

  from OP_waveforms import *

  path = '../data/Ijen'
  list_sta = ['DAM','IJEN','KWUI','MLLR','POS','POSI','PSG','RAUN','TRWI']
  for sta in list_sta:
    list_files = glob.glob(os.path.join(path,'*%s*Z*'%sta))
    list_files.sort()
    for file in list_files:
      if os.path.isfile(file):
        print file
        wf = Waveform()
        wf.read_from_file(file)
        wf.rmean()
        wf.bp_filter(1,10)
        x = wf.values

        #kurt_file = '%s/kurt/%s_KURT'%(os.path.dirname(file),os.path.basename(file))
        kurt_file = '/media/disk/location/%s_KURT'%os.path.basename(file)
        if not os.path.isfile(kurt_file):
          print "computing kurtosis..."
          kurt_window = 3.
          wf.process_kurtosis(kurt_window,recursive=True,pre_taper=True,post_taper=True)
          #wf.values[:1700] = 0
          wf.write_to_file_filled(kurt_file,format='MSEED',fill_value=0)
          x_kurt = wf.values
          os.symlink(kurt_file,os.path.join(path,os.path.basename(kurt_file)))

        #grad_file = '%s/grad/%s_GRAD'%(os.path.dirname(file),os.path.basename(file))
        grad_file = '/media/disk/location/%s_GRAD'%os.path.basename(file)
        if not os.path.isfile(grad_file):
          print "computing first derivative"
          wf.take_positive_derivative(pre_taper=True,post_taper=True)
          wf.write_to_file_filled(grad_file,format='MSEED',fill_value=0)
          x_grad = wf.values
          os.symlink(grad_file,os.path.join(path,os.path.basename(grad_file)))

        #plotting_traces(x,x_kurt,x_grad)



def plotting_traces(x,y,z):

  import matplotlib.pyplot as plt
  fig = plt.figure()
  fig.set_facecolor('white')
  ax1 = fig.add_subplot(311)
  ax1.plot(x,'k')
  ax2 = fig.add_subplot(312)
  ax2.plot(y,'k')
  ax3 = fig.add_subplot(313)
  ax3.plot(z,'k')
  plt.show()


def sds_link():

  data_path = '/media/disk/location'
  dest_path = '../data/Ijen'
  list_sta = ['DAM','IJEN','KWUI','MLLR','POS','POSI','PSG','RAUN','TRWI']
  for sta in list_sta:
    station = sta.lower()
    list_files = glob.glob(os.path.join(data_path,station,'*Z*'))
    list_files.sort()
    for file in list_files:
      dest_file = os.path.join(dest_path,os.path.basename(file))
      os.symlink(file,dest_file)


def migration_and_location(opdict):

  from OP_waveforms import *
  from NllGridLib import read_stations_file, read_hdr_file
  import h5py
  from migration import migrate_4D_stack, extract_max_values
  from locations_trigger import trigger_locations_inner, number_good_kurtosis_for_location,write_header_options

  base_path = opdict['base_path']
  data_dir = os.path.join(base_path,'data',opdict['datadir'])

  dataglob = opdict['dataglob']

  grid_filename_base = os.path.join(base_path,'lib',opdict['time_grid'])
  search_grid_filename = os.path.join(base_path,'lib',opdict['search_grid'])
  grid_info = read_hdr_file(search_grid_filename)

  nx = grid_info['nx']
  ny = grid_info['ny']
  nz = grid_info['nz']
  n_buf = nx*ny*nz

  npts = int(opdict['data_length']*opdict['fs'])
  delta = 1./opdict['fs']
  loclevel = opdict['loclevel']

  snr_limit = opdict['snr_limit']
  snr_tr_limit = opdict['snr_tr_limit']
  sn_time = opdict['sn_time']
  n_kurt_min = opdict['n_kurt_min']

  from hdf5_grids import get_interpolated_time_grids
  time_grids = get_interpolated_time_grids(opdict)

  import pandas as pd
  df = pd.read_csv('../lib/Ijen_list_vt.csv',index_col=False)
  dates = df.index

  n_ok = 0
  loc_filename = os.path.join(base_path,'out',opdict['outdir'],'loc','locations.dat')
  loc_file = open(loc_filename,'w')
  write_header_options(loc_file,opdict)

  for date in dates:
    data = {}
    data_files = glob.glob(os.path.join(data_dir,'*%s_%s*%s*'%(str(date)[:8],str(date)[8:],opdict['dataglob'])))
    kurt_files = glob.glob(os.path.join(data_dir,'*%s_%s*%s*'%(str(date)[:8],str(date)[8:],opdict['kurtglob'])))
    grad_files = glob.glob(os.path.join(data_dir,'*%s_%s*%s*'%(str(date)[:8],str(date)[8:],opdict['gradglob'])))
    data_files.sort()
    kurt_files.sort()
    grad_files.sort()
    filename = '%s_%d.hdf5'%(opdict['outdir'],date)
    grid_file = os.path.join(base_path,'out',opdict['outdir'],'grid',filename)
    stack_file = os.path.join(base_path,'out',opdict['outdir'],'stack','stack_all_'+filename)
    for file in grad_files:
      wf = Waveform()              
      wf.read_from_file(file)
      data[wf.station] = wf.values
      start_time = wf.starttime

    logging.info('Doing migration to %s'%grid_file)
    f = h5py.File(grid_file,'w')
    stack_grid = f.create_dataset('stack_grid',(n_buf,npts),'f',chunks=(1,npts))
    stack_shift_time = migrate_4D_stack(data,delta,time_grids,stack_grid)
    n_buf,nt = stack_grid.shape

    # add useful information to dataset
    for key,value in grid_info.iteritems():
      stack_grid.attrs[key] = value
    stack_grid.attrs['dt'] = delta
    stack_grid.attrs['start_time'] = -stack_shift_time

    # extract max-stack
    logging.info('Extracting max_val etc. to %s'%stack_file)
    f_stack = h5py.File(stack_file,'w')
    # extract maxima
    extract_max_values(stack_grid,grid_info,f_stack)
    for name in f_stack:
      dset = f_stack[name]
      logging.debug('After extract_max_values : %s %f %f'%(name,np.max(dset),np.sum(dset)))
      dset.attrs['start_time'] = -stack_shift_time
      dset.attrs['dt'] = delta

    max_val = f_stack['max_val']
    max_x = f_stack['max_x']
    max_y = f_stack['max_y']
    max_z = f_stack['max_z']
    max_val_smoothed = smooth(max_val[:])

    locs = trigger_locations_inner(max_val_smoothed,max_x,max_y,max_z,loclevel,loclevel,-stack_shift_time,delta)

    # close the stack and grid files 
    f_stack.close()
    f.close()
    logging.info('Saved 4D grid to file %s'%grid_file)

    for loc in locs:
      loc['o_time'] = start_time + loc['o_time']
      if number_good_kurtosis_for_location(kurt_files,data_files,loc,time_grids,snr_limit,snr_tr_limit,sn_time) > n_kurt_min:
        logging.info("Max = %.2f, %s - %.2fs + %.2f s, x=%.4f pm %.4f km, y=%.4f pm %.4f km, z=%.4f pm %.4f km"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
        loc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
        n_ok=n_ok+1

  loc_file.close()



def location_only(opdict):

  from locations_trigger import trigger_locations_inner, number_good_kurtosis_for_location,write_header_options
  from h5py

  base_path = opdict['base_path']

  snr_limit = opdict['snr_limit']
  snr_tr_limit = opdict['snr_tr_limit']
  sn_time = opdict['sn_time']
  n_kurt_min = opdict['n_kurt_min']

  dt = 1./optdict['fs']

  n_ok = 0
  loc_filename = os.path.join(base_path,'out',opdict['outdir'],'loc','locations.dat')
  loc_file = open(loc_filename,'w')
  write_header_options(loc_file,opdict)

  stack_files = glog.glob(os.path.join(base_path,'out',opdict['outdir'],'stack','*'))
  stack_files.sort()
  for stack_filename in stack_files:
    f_stack = h5py.File(stack_filename,'r')
    max_val = f_stack['max_val']
    max_x = f_stack['max_x']
    max_y = f_stack['max_y']
    max_z = f_stack['max_z']
    max_val_smoothed = smooth(max_val[:])

    locs = trigger_locations_inner(max_val_smoothed,max_x,max_y,max_z,loclevel,loclevel,max_val.attrs['start_time'],dt)

    f_stack.close()

    for loc in locs:
      loc['o_time'] = start_time + loc['o_time']
      if number_good_kurtosis_for_location(kurt_files,data_files,loc,time_grids,snr_limit,snr_tr_limit,sn_time) > n_kurt_min:
        logging.info("Max = %.2f, %s - %.2fs + %.2f s, x=%.4f pm %.4f km, y=%.4f pm %.4f km, z=%.4f pm %.4f km"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
        loc_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(loc['max_trig'],loc['o_time'].isoformat(),loc['o_err_left'], loc['o_err_right'],loc['x_mean'],loc['x_sigma'],loc['y_mean'],loc['y_sigma'],loc['z_mean'],loc['z_sigma']))
        n_ok=n_ok+1

  loc_file.close()



if __name__ == '__main__':

  #list_of_files()
  #ijen()
  #sds_link()

  from options import WavelocOptions
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')

  wo = WavelocOptions()
  args = wo.p.parse_args()

  wo.set_options()
  wo.verify_migration_options()
  wo.verify_location_options()

  migration_and_location(wo.opdict)
  location_only(wo.opdict)
