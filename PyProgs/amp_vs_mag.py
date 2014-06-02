#!/usr/bin/env python

import os,glob,sys
import numpy as np
from locations_trigger import read_locs_from_file
import matplotlib.pyplot as plt


def amp_vs_mag(crises):

  amp_all, mag_all = [],[]

  for date in crises:

    print date
    loc_filename="../out/%s/loc/locations.dat"%date
    locs=read_locs_from_file(loc_filename)

    for loc in locs:
      mag_all.append(loc['ml'])
      amp_all.append(loc['max_trig'])


  fig = plt.figure()
  fig.set_facecolor('white')
  plt.plot(amp_all,mag_all,'kx')
  plt.xlabel('Maximum stack amplitude')
  plt.ylabel('Magnitude')
  plt.savefig('/home/nadege/Desktop/amp_vs_mag.pdf')
  plt.show()


def time_vs_mag(crises):

  import matplotlib.gridspec as gridspec

  id = 0
  fig = plt.figure(figsize=(15,8))
  fig.set_facecolor('white')
  G = gridspec.GridSpec(6,2)
  type=['i','i','i','e','e','i','e','i','e','i+e','i']

  for date in crises:

    otime_all, mag_all = [],[]

    print date
    id += 1
    loc_filename="../out/%s/loc/locations.dat"%date
    locs=read_locs_from_file(loc_filename)

    for loc in locs:
      mag_all.append(loc['ml'])
      otime_all.append(float(loc['o_time']))

    if id <= 6:
      ax = fig.add_subplot(G[id-1,0])
    else:
      ax = fig.add_subplot(G[id-7,1])
    #ax.set_axis_off()
    ax.plot((otime_all-np.min(otime_all))/1000,mag_all,'k*',markersize=1)
    ax.text(0.78,0.78,"%s %s"%(type[id-1],date),transform=ax.transAxes)

    if id == len(crises):
      ax.set_xlabel('Time sample/1000')

  plt.savefig('../../Desktop/mag_vs_time.pdf')
  plt.show()


if __name__ == '__main__':

  crises=['2009-10-14','2009-10-18','2009-10-29','2009-11-05','2009-12-14','2009-12-29','2010-01-02','2010-09-23','2010-10-14','2010-12-09','2011-02-02']
  #amp_vs_mag(crises)
  time_vs_mag(crises)


