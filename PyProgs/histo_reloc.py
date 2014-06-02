#!/usr/bin/env python
# encoding: utf-8

import glob,os,sys
import numpy as np
import matplotlib.pyplot as plt
from NllGridLib import qd_read_hyp_file, qd_read_picks_from_hyp_file
from obspy.core import utcdatetime
import matplotlib.mlab as mlab
from scipy.interpolate import UnivariateSpline
from locations_trigger import read_locs_from_file, read_header_from_file

def histo(locdir,loca):

  common_file = open('../out/%s/common_reloc.dat'%locdir,'w')

  # Manual locations
  hyp_path="../data/14_oct_locpicks/loc_hyp"
  hyp_files=glob.glob("%s/piton*.hyp"%hyp_path)
  hyp_parameters=[]
  recap=open('../data/14_oct_locpicks/loc_hyp/recap','w')
  for hyp_file in hyp_files:
    (otime,hypo_x,sigma_x,hypo_y,sigma_y,hypo_z,sigma_z)=qd_read_hyp_file(hyp_file)
    phase_dict=qd_read_picks_from_hyp_file(hyp_file)
    hyp_parameters.append(((otime,hypo_x,sigma_x,hypo_y,sigma_y,hypo_z,sigma_z),phase_dict))
    recap.write('%s, x= %s km, y= %s km, z= %s km\n'%(otime, hypo_x, hypo_y, hypo_z))
  number_hyp=len(hyp_files)
  recap.close()

  dx,dy,dz,dt,d,inc=[],[],[],[],[],[]
  common = 0

  locpath="%s/out/%s"%(path,locdir)
  loc_filename="%s/%s"%(locpath,loca)

  # Read location file
  locs=read_locs_from_file(loc_filename)
  head=read_header_from_file(loc_filename,{})
  number_of_locations=len(locs)

  for loc in locs:
    stack_max=loc['max_trig']
    stack_time=loc['o_time']
    stack_time_err_left=loc['o_err_left']
    stack_time_err_right=loc['o_err_right']
    stack_x=loc['x_mean']
    stack_x_err=loc['x_sigma']
    stack_y=loc['y_mean']
    stack_y_err=loc['y_sigma']
    stack_z=loc['z_mean']
    stack_z_err=loc['z_sigma']

    A=[(np.abs(hyp_tuple[0][0]-stack_time),hyp_tuple) for hyp_tuple in hyp_parameters ]
    A.sort()
    hyp_tuple=A[0][1]
    (otime, hypo_x, sigma_x, hypo_y, sigma_y, hypo_z, sigma_z)=hyp_tuple[0]
    phases=hyp_tuple[1]

    if np.abs(otime-stack_time) < 1.0:
      common+=1
      common_file.write('%s, x= %s km, y= %s km, z= %s km\n'%(otime, hypo_x, hypo_y, hypo_z))
      dt.append(stack_time-otime)
      dz.append(-stack_z+hypo_z)
      dy.append(stack_y-hypo_y)
      dx.append(stack_x-hypo_x)
      dauto=np.sqrt(stack_x**2+stack_y**2+stack_z**2)
      dman=np.sqrt(hypo_x**2+hypo_y**2+hypo_z**2)
      dist=np.sqrt((stack_x-hypo_x)**2+(stack_y-hypo_y)**2+(-stack_z+hypo_z)**2)
      d.append(dist)
      inc.append(1/np.sqrt(stack_x**2+stack_y**2+stack_z**2)*np.sqrt((stack_x**2*stack_x_err**2+stack_y**2*stack_y_err**2+stack_z**2*stack_y_err**2)))

  common_file.close()

  comment=''
  if head['krec'] == 'True':
    comment='rec'
  if head['kderiv'] == 'True':
    comment='%s deriv'%comment
  if head['gauss'] == 'True':
    comment='%s gauss'%comment
  if head['reloc'] == 'True':
    comment='%s reloc'%comment

  spec = r'%.0f-%.0fHz %s ($K_t$=%s)'%(head['c1'],head['c2'],comment,head['loclevel'])

  mt = np.mean(dt)
  mx = np.mean(dx)
  my = np.mean(dy)
  mz = np.mean(dz)

  lab={}
  lab['time']='%s \n(W: %d - C: %d - M: %.2f)'%(spec,number_of_locations,common,mt)
  lab['x']='%s \n(W: %d - C: %d - M: %.2f)'%(spec,number_of_locations,common,mx)
  lab['y']='%s \n(W: %d - C: %d - M: %.2f)'%(spec,number_of_locations,common,my)
  lab['z']='%s \n(W: %d - C: %d - M: %.2f)'%(spec,number_of_locations,common,mz)

  plot(number_of_locations,common,number_hyp,mt,dt,lab,'time')
  plot(number_of_locations,common,number_hyp,mx,dx,lab,'x')
  plot(number_of_locations,common,number_hyp,my,dy,lab,'y')
  plot(number_of_locations,common,number_hyp,mz,dz,lab,'z')

# -------------------------------------------------------------------------- #
def spline(bins,n):
  c='k'
  linestyle='-'
  lw=3.5
  s1=UnivariateSpline(bins[0:-1],n)
  s2=UnivariateSpline(bins,s1(bins))
  y=s2(bins)
  maxi=np.argmax(y)
  plt.plot(bins,y,color=c,lw=lw,ls=linestyle)

# -------------------------------------------------------------------------- #
def plot(number_of_locations,common,number_hyp,m,d,lab,comp):

  c = 'k'

  if comp == 'x':
    let='(b)'
  elif comp == 'y':
    let='(c)'
  else:
    let='(d)'

  fig=plt.figure(figsize=(11,9))
  fig.set_facecolor('white')
  if comp ==  'time':
    int=[-1,1]
    nb = 50
    n, bins, patches = plt.hist(d,nb,range=None,align='left',histtype='stepfilled',alpha=.2,color=c,label=lab[comp])#,fill=False,lw=0)
    plt.figtext(0.05,0.92,"(a)",fontsize=24)
    plt.legend(loc=2,prop={'size':'x-large'})
    plt.axis(int+[0,np.max(n)+3])
    spline(bins,n)
    unit='s'
    plt.xticks(size='x-large')
    plt.yticks(size='x-large')

  else:
    int=[-3,3]
    nb = 100
    n, bins, patches = plt.hist(d,nb,range=None,align='left',histtype='stepfilled',alpha=.2,color=c,label=lab[comp])#,fill=False,lw=0)
    plt.figtext(0.05,0.92,let,fontsize=24)
    plt.axis(int+[0,np.max(n)+3])
    plt.legend(loc=2,prop={'size':'x-large'})
    spline(bins,n)
    unit='km'
    plt.xticks(size='x-large')
    plt.yticks(size='x-large')

  plt.xlabel('$\Delta$%s (%s)'%(comp,unit),fontsize=24)
  plt.ylabel('Number of common events',fontsize=24)
  plt.title('$\Delta$%s'%comp,fontsize=24)

#########################################################################
locdirs="2010-10-14/loc"
histo(locdirs,'relocations.dat')

#plt.show()
