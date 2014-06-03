import os, glob, unittest, h5py
import numpy as np
from options import WavelocOptions
from OP_waveforms import Waveform
from migration import do_migration_setup_and_run
from integrate4D import * 
from synth_migration import generateSyntheticDirac
from test_processing import waveforms_to_signature



def suite():
  suite = unittest.TestSuite()
  suite.addTest(SyntheticMigrationTests('test_dirac_migration'))
  #suite.addTest(MigrationTests('test_migration'))
  #suite.addTest(MigrationTests('test_migration_fullRes'))
  return suite

def hdf5_to_signature(base_path,datadir,dataglob,output_filename):

  sig_file=open(os.path.join(base_path,datadir,output_filename),'w')
  allfiles=glob.glob(os.path.join(base_path,datadir, dataglob))
  for filename in allfiles :
    basename=os.path.basename(filename)
    f=h5py.File(filename,'r')
    for name in f:
      logging.debug('Signature for %s %s : '%(basename,name))
      dset=f[name]
      maximum=np.max(dset)
      datasum=np.sum(dset)
      sig_file.write("%s \t %s \t %.6f \t %.6f\n"%(basename,name,maximum,datasum))
    f.close()


class SyntheticMigrationTests(unittest.TestCase):

  def test_dirac_migration(self):
    from locations_trigger import trigger_locations_inner
    from filters import smooth

    wo=WavelocOptions()
    wo.set_test_options()
    
    wo.opdict['outdir'] = 'TEST_Dirac'
    wo.opdict['load_ttimes_buf'] = True # Optimized in time, but you must be sure you're reading the right grid for the test
    wo.opdict['syn_addnoise'] = False

    # Synthetics for IJEN network
    wo.opdict['syn_snr'] = 10
    wo.opdict['syn_amplitude'] = [5.0]
    wo.opdict['syn_samplefreq'] = 100.0
    wo.opdict['syn_kwidth'] = 0.1
    wo.opdict['syn_otime'] = [10.0]
    wo.opdict['loclevel'] = 3
    wo.opdict['stations'] = '/home/nadege/waveloc/lib/coord_stations_ijen'
    wo.opdict['sta_list'] = "IJEN,TRWI,RAUN,MLLR,KWUI,POSI,POS,PSG,DAM"
    wo.opdict['search_grid'] = 'grid.ijen.search.hdr'
    wo.opdict['syn_filename'] = 'Ijen_test_grid4D_hires.hdf5'
    wo.opdict['time_grid'] = 'ijen.P'
    wo.opdict['syn_datalength']= 40.0

    nbx,nby,nbz = 41,41,11
    nbuf = nbx*nby*nbz

    plot = True
    save = False
    if save:
      f_unc = h5py.File('/home/nadege/waveloc/out/TEST_Dirac/incertitude.hdf5','w')
      unc_grid_x = f_unc.create_dataset('unc_grid_x',(nbuf,))
      unc_grid_y = f_unc.create_dataset('unc_grid_y',(nbuf,))
      unc_grid_z = f_unc.create_dataset('unc_grid_z',(nbuf,))

      #err_grid_x = f_unc.create_dataset('err_grid_x',(nbuf,))
      #err_grid_y = f_unc.create_dataset('err_grid_y',(nbuf,))
      #err_grid_z = f_unc.create_dataset('err_grid_z',(nbuf,))

    if plot:
      # uncertainty equal to 0
      #list_ibb = [382,404,1316,1767,2207,2658,3098,4890,5307,5803,5867,5868,5924,6376,7225,7276,7288,7673,7674,7675,7684,7685,7696,7728,7729,7739,7740,7741,7751,8123,8124,8125,8126,8127,8134,8145,8146,8147,8149,8157,8168,8212,8573,8574,8575,8576,8577,8578,8585,8587,8608,8618,8619,9024,9025,9026,9027,9028,9029,9030,9037,9047,9048,9476,9477,9480,9481,9488,9489,9490,9491,9498,9499,9521,9532,9544,9554,9928,9929,9930,9931,9932,9934,9938,9940,9941,9952,9972,9973,9983,10005,10379,10380,10381,10382,10383,10385,10386,10403,10411,10412,10423,10424,10829,10830,10831,10832,10837,10851,10852,10862,10875,10885,11282,11283,11297,11304,11557,11732,11733,11734,11743,11755,11776,11777,11997,12183,12184,12195,12205,12635,12765,12904,18172]

      # maximum uncertainty values
      list_ibb = [1794 ,9716,11681,12234,10178,12683,16462,13167,11738,16809,18240,12091,10306,18304,13258,15845,16381,9933,9727,16830,13123,16340,9264,12004,15070,15346,11609,12544,11793,17327,17732,13553,18260,12858,18250,17359,10168,9717,12190,9255,10619,16436,15556,10835,17326,12287,14643,16477,10640,10639,17336,15961,11091,12255,16866,18305,12531,12543,13652,17788,12688,17270,15895,16775,12247,12475,17778,6182,10179,16380,18249,15951,15564,18239,10630,12188,15488,11726,12476,11147]

    for ixx in range(0,nbx):
      #if ixx != 15:
      #  listiyy = [15]
      #else:
      #  listiyy = range(0,31)
      #for iyy in listiyy:
      for iyy in range(0,nby):
        for izz in range(0,nbz):
          wo.opdict['syn_ix'] = [ixx]
          wo.opdict['syn_iy'] = [iyy]
          wo.opdict['syn_iz'] = [izz]

          ibb = ixx*nby*nbz + iyy*nbz + izz
          if plot:
            if ibb not in list_ibb:
              continue
          print "\nINDICE : ",ibb, ixx, iyy, izz

          wo.verify_migration_options()
          wo.verify_location_options()
          wo.verify_synthetic_options()

    ##########################
    # generate the test case and retrieve necessary information
    ##########################

          logging.info('Running synthetic test case generation...')
          test_info = generateSyntheticDirac(wo.opdict)
          logging.debug(test_info)

          # retrieve info
          grid_filename = test_info['dat_file']
          stack_filename = test_info['stack_file']
          nx,ny,nz,nt = test_info['grid_shape']
          dx,dy,dz,dt = test_info['grid_spacing']
          x_orig,y_orig,z_orig = test_info['grid_orig']
          ix_true,iy_true,iz_true,it_true = test_info['true_indexes']
          stack_start_time = test_info['start_time']

          # plot base filename
          base_path = wo.opdict['base_path']
          outdir = wo.opdict['outdir']

          # loclevel for triggers
          loclevel = wo.opdict['loclevel']

          # set up x, y, z, t arrays
          x = np.arange(nx)*dx
          y = np.arange(ny)*dy
          z = np.arange(nz)*dz
          t = np.arange(nt)*dt+stack_start_time

          # extract the max stacks
          f_stack = h5py.File(stack_filename,'r')
          max_val = f_stack['max_val_smooth']
          max_x = f_stack['max_x']
          max_y = f_stack['max_y']
          max_z = f_stack['max_z']

          locs = trigger_locations_inner(max_val,max_x,max_y,max_z,loclevel,loclevel,stack_start_time,dt)
          print locs

          x_found = round(locs[0]['x_mean'],1)
          y_found = round(locs[0]['y_mean'],1)
          z_found = round(locs[0]['z_mean'],1)
          print wo.opdict['syn_otime'],locs[0]['o_time']
          print np.array(wo.opdict['syn_ix'])*dx+x_orig,x_found
          print np.array(wo.opdict['syn_iy'])*dy+y_orig,y_found
          print np.array(wo.opdict['syn_iz'])*dz+z_orig,z_found

          print "UNCERTAINTIES :",locs[0]['x_sigma'],locs[0]['y_sigma'],locs[0]['z_sigma']
          test_info['o_time'] = locs[0]['o_time']
          test_info['x_err'] = (locs[0]['x_mean']-locs[0]['x_sigma'],locs[0]['x_mean']+locs[0]['x_sigma'])
          test_info['y_err'] = (locs[0]['y_mean']-locs[0]['y_sigma'],locs[0]['y_mean']+locs[0]['y_sigma'])
          test_info['z_err'] = (locs[0]['z_mean']-locs[0]['z_sigma'],locs[0]['z_mean']+locs[0]['z_sigma'])
          test_info['t_err'] = (locs[0]['o_time']-locs[0]['o_err_left'],locs[0]['o_time']+locs[0]['o_err_right'])
          #test_info['t_err'] = (locs[0]['o_time']-locs[0]['o_err_left']-stack_start_time,locs[0]['o_time']+locs[0]['o_err_right']-stack_start_time)

          if plot:

            # The first figure shows the stacking values for the true origin time.
            # The second figure shows the stacking values for the recovered origin time.
            # The third plot shows the kurtosis traces alignment.
            # yellow star = true location ; red star = Waveloc's location

            it_found = wo.opdict['syn_samplefreq'] * test_info['o_time']
            ix_found = int(round((x_found-x_orig)*1./dx))
            iy_found = int(round((y_found-y_orig)*1./dy))
            iz_found = int(round((z_found-z_orig)*1./dz))

            from plot_mpl import plotDiracTest
            print test_info['true_indexes']
            plotDiracTest(test_info,'/home/nadege/waveloc/out/TEST_Dirac/fig',2,ixx,iyy,izz,ix_found,iy_found,iz_found)

            test_info['true_indexes'] = (ix_found, iy_found, iz_found, it_found)
            print test_info['true_indexes']
            plotDiracTest(test_info,'/home/nadege/waveloc/out/TEST_Dirac/fig',2,ixx,iyy,izz,ix_found,iy_found,iz_found)

            wix = np.argmin(np.abs(np.array(wo.opdict['syn_ix'])*dx+x_orig - np.arange(x_orig,x_orig+nx*dx,dx)))
            wiy = np.argmin(np.abs(np.array(wo.opdict['syn_iy'])*dy+y_orig - np.arange(y_orig,y_orig+ny*dy,dy)))
            wiz = np.argmin(np.abs(np.array(wo.opdict['syn_iz'])*dz+z_orig - np.arange(z_orig,z_orig+nz*dz,dz)))
            wib = wix*ny*nz + wiy*nz + wiz

            from hdf5_grids import get_interpolated_time_grids
            time_grids = get_interpolated_time_grids(wo.opdict)
            ttimes = {}
            for sta in wo.opdict['sta_list'].split(','):
              if time_grids.has_key(sta):
                ttimes[sta] = time_grids[sta].grid_data[wib]
              else:
                logging.info('Missing travel-time information for station %s. Ignoring station...'%sta)

            import matplotlib.pyplot as plt
            fig = plt.figure()
            fig.set_facecolor('white')
            n_traces = len(ttimes) + 2
            list_ind = [] 
            for iii,key in enumerate(sorted(ttimes)):
              delta_t = test_info['grid_spacing'][-1]
              t = np.arange(0,len(test_info['data'][key])*delta_t,delta_t)
              ax = fig.add_subplot(n_traces,2,2*iii+1)
              ax.set_axis_off()
              ax.plot(t,test_info['data'][key],'k')
              ax.text(0.2,0.5,key)

              ind = np.argmin(np.abs(t-ttimes[key]))
              t = np.arange(0,len(test_info['data'][key][ind:ind+len(max_val)])*delta_t,delta_t)
              ax = fig.add_subplot(n_traces,2,2*iii+2)
              ax.set_axis_off()
              if len(t) == len(test_info['data'][key][ind:ind+len(max_val)]):
                ax.plot(t,test_info['data'][key][ind:ind+len(max_val)],'k')
              else:
                ax.plot(t[:-1],test_info['data'][key][ind:ind+len(max_val)],'k')

            t = np.arange(0,len(max_val)*delta_t,delta_t)
            ind = np.argmin(np.abs(t-locs[0]['o_time']))
            ax = fig.add_subplot(n_traces,2,2*(iii+1)+2)
            ax.set_axis_off()
            t = np.arange(0,len(max_val[ind:])*delta_t,delta_t)
            ax.plot(max_val,'k')
            #if len(t) == len(max_val[ind:]):
            #  ax.plot(t,max_val[ind:],'k')
            #else:
            #  ax.plot(t[:-1],max_val[ind:],'k')
            plt.show()

          if save:
            #f_unc['err_grid_x'][ibb] = np.abs(locs[0]['x_mean']-(wo.opdict['syn_ix'][0]*dx+x_orig))
            #f_unc['err_grid_y'][ibb] = np.abs(locs[0]['y_mean']-(wo.opdict['syn_iy'][0]*dy+y_orig))
            #f_unc['err_grid_z'][ibb] = np.abs(locs[0]['z_mean']-(wo.opdict['syn_iz'][0]*dz+z_orig))
            f_unc['unc_grid_x'][ibb] = locs[0]['x_sigma']
            f_unc['unc_grid_y'][ibb] = locs[0]['y_sigma']
            f_unc['unc_grid_z'][ibb] = locs[0]['z_sigma']

          f_stack.close()

    if save:
      f_unc.close()


class MigrationTests(unittest.TestCase):

  def setUp(self):

    self.wo=WavelocOptions()
    self.wo.set_test_options()
    self.wo.verify_migration_options()



#  @unittest.skip('Not running small test')
  @unittest.expectedFailure
  def test_migration(self):

    self.wo.opdict['load_ttimes_buf'] = True
    self.wo.opdict['data_length'] = 300

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    expected_signature_filename = os.path.join(base_path,test_datadir,'test_stack_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()


    do_migration_setup_and_run(self.wo.opdict)

    #waveforms_to_signature(base_path,os.path.join('out',outdir,'stack'),'stack*mseed','stack_signature.dat')
    hdf5_to_signature(base_path,os.path.join('out',outdir,'stack'),'stack_all_2010-10-14T00:14:00.000000Z.hdf5','stack_signature.dat')
    signature_filename=os.path.join(base_path,'out',outdir,'stack','stack_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)

  @unittest.skip('Not running full resolution test')
  #@profile
  def test_migration_fullRes(self):

    self.wo.opdict['search_grid'] = 'grid.Taisne.search.hdr'
    self.wo.opdict['outdir'] = 'TEST_fullRes'
    self.wo.opdict['load_ttimes_buf'] = True
    self.wo.opdict['data_length'] = 300
    self.wo.verify_migration_options()

    base_path=self.wo.opdict['base_path']
    test_datadir=self.wo.opdict['test_datadir']
    outdir=self.wo.opdict['outdir']

    expected_signature_filename = os.path.join(base_path,test_datadir,'TEST_fullRes_stack_signature.dat')
    expected_signature_file = open(expected_signature_filename,'r') 
    expected_lines=expected_signature_file.readlines()

    do_migration_setup_and_run(self.wo.opdict)

    #waveforms_to_signature(base_path,os.path.join('out',outdir,'stack'),'stack*mseed','stack_signature.dat')
    hdf5_to_signature(base_path,os.path.join('out',outdir,'stack'),'stack*hdf5','stack_signature.dat')
    signature_filename=os.path.join(base_path,'out',outdir,'stack','stack_signature.dat')
    signature_file = open(signature_filename,'r') 
    lines=signature_file.readlines()

    self.assertSequenceEqual(lines,expected_lines)



if __name__ == '__main__':

  import logging
  logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(asctime)s : %(message)s')
 
  unittest.TextTestRunner(verbosity=2).run(suite())
 
