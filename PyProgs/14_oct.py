#!/usr/bin/env python
# encoding: utf-8

import os, sys, optparse, glob
from OP_waveforms import *
from obspy.core import utcdatetime
import time


def piton_de_la_fournaise():

  datadir='/home/nadege/waveloc/data/test_recurs'
  #datadir='/home/nadege/Desktop/raw_data_2_crises/309'

  #dataglob="*"
  #filtglob='*filt.mseed'

  kurtglob="*filt_kurt3.sac"
  dataglob=kurtglob

  gradglob="*filt_kurt3_grad.sac"
  gaussglob="*filt_kurt3_grad_gauss.sac"

  kurt_window=1

  datafiles=glob.glob(os.path.join(datadir,dataglob))
  datafiles.sort()

  t1=utcdatetime.UTCDateTime("2010-10-14T00:00:00.0Z")
  t2=utcdatetime.UTCDateTime("2010-10-14T16:00:00.0Z")

  for datafile in datafiles:
    wf=Waveform()
    wf.read_from_file(datafile,starttime=t1,endtime=t2)
    sta=wf.station
    cha=wf.channel

    print sta

    #print "filtering..."
    #if sta != 'UV01':
      #wf.bp_filter(20,35,rmean=True,taper=True)
    #else:
      #wf.bp_filter(4,10,rmean=True,taper=True)
    #wf.write_to_file_filled('/home/nadege/waveloc/data/Piton/2009-11-05/%s%s'%(os.path.basename(datafile.split(dataglob[1:])[0]),filtglob[1:]),format='SAC',fill_value=0)
    #wf.write_to_file_filled('/home/nadege/waveloc/data/Piton/2009-11-05/%s/2009-11-05T00:00:00.YA.%s.%s.filt.mseed'%(cha,sta,cha),format='MSEED',fill_value=0)

    #tref=time.time()
    #print "computing kurtosis..."
    #wf.process_kurtosis(kurt_window,recursive=True,pre_taper=True,post_taper=True)
    #wf.write_to_file_filled('/home/nadege/waveloc/data/Piton/2009-11-05/s%s'%(os.path.basename(datafile.split(dataglob[1:])[0]),kurtglob[1:]),format='SAC',fill_value=0)
    #wf.write_to_file_filled('/home/nadege/waveloc/data/Piton/2009-11-05/%s/2009-11-05T00:00:00.YA.%s.%s.filt_kurt.mseed'%(cha,sta,cha),format='MSEED',fill_value=0)
    #print time.time()-tref

    print "computing first derivative"
    wf.take_positive_derivative(pre_taper=True,post_taper=True)
    wf.write_to_file_filled('/home/nadege/waveloc/data/test_recurs/%s%s'%(os.path.basename(datafile.split(dataglob[1:])[0]),gradglob[1:]),format='SAC',fill_value=0)
    #wf.write_to_file_filled('/home/nadege/waveloc/data/Piton/2009-11-05/%s/2009-11-05T00:00:00.YA.%s.%s.filt_kurt_grad.mseed'%(cha,sta,cha),format='MSEED',fill_value=0)

    print "replacing kurtosis by gaussian distributions"
    wf.process_gaussian(10,0,.1)
    wf.write_to_file_filled('/home/nadege/waveloc/data/test_recurs/%s%s'%(os.path.basename(datafile.split(dataglob[1:])[0]),gaussglob[1:]),format='SAC',fill_value=0)
    #wf.write_to_file_filled('/home/nadege/waveloc/data/Piton/2009-11-05/%s/2009-11-05T00:00:00.YA.%s.%s.filt_kurt_grad_gauss.mseed'%(cha,sta,cha),format='MSEED',fill_value=0)


def ijen():

  #path = '/media/disk/CLASSIFICATION'
  path = '/media/disk'
  #list_sta = glob.glob(os.path.join(path,'*'))
  #list_sta.sort()
  list_sta = ['%s/psg'%path]
  for sta_path in list_sta:
    list_files = glob.glob(os.path.join(sta_path,'ID.*'))
    list_files.sort()
    for file in list_files:
      filt_file = '%s/%s_FILT'%(os.path.dirname(file),os.path.basename(file))
      kurt_file = '%s/kurt/%s_KURT'%(os.path.dirname(file),os.path.basename(file))
      grad_file = '%s/grad/%s_GRAD'%(os.path.dirname(file),os.path.basename(file))
      if os.path.isfile(file):
        print file
        wf = Waveform()
        wf.read_from_file(file)
        wf.rmean()
        wf.bp_filter(1,10)
        x = wf.values
        #wf.write_to_file_filled(filt_file,format='MSEED',fill_value=0)

        #kurt_file = '/home/nadege/Desktop/kurt/%s_KURT'%os.path.basename(file)
        if not os.path.isfile(kurt_file):
          print "computing kurtosis..."
          kurt_window = 3.
          wf.process_kurtosis(kurt_window,recursive=True,pre_taper=True,post_taper=True)
          #wf.values[:1700] = 0
          wf.write_to_file_filled(kurt_file,format='MSEED',fill_value=0)
          x_kurt = wf.values

        #grad_file = '/home/nadege/Desktop/grad/%s_GRAD'%os.path.basename(file)
        if not os.path.isfile(grad_file):
          print "computing first derivative"
          wf.take_positive_derivative(pre_taper=True,post_taper=True)
          wf.write_to_file_filled(grad_file,format='MSEED',fill_value=0)
          x_grad = wf.values

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


if __name__ == '__main__':
  #piton_de_la_fournaise()
  ijen()

