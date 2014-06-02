#!/usr/bin/env python

import os, glob,sys
import numpy as np
import time
from scipy.io.matlab import mio
from mayavi import mlab
from obspy.core import utcdatetime
import matplotlib.pyplot as plt
from matplotlib import mpl
import math
# -----------------------------------------------------------------------
def plot_map_3d(matfile):
  mat=mio.loadmat(matfile)
  XIsub=mat["XIsub"][0]/1000.0
  YIsub=mat["YIsub"][:,0]/1000.0
  ZIsub=np.matrix(mat["ZIsub"]/1000.0)

  xmin=np.min(XIsub)
  xmax=np.max(XIsub)
  ymin=np.min(YIsub)
  ymax=np.max(YIsub)
  zmin=np.min(ZIsub)
  zmax=np.max(ZIsub)
  ext=[xmin,xmax,ymin,ymax,zmin,zmax]

  return XIsub,YIsub,ZIsub,ext
# -----------------------------------------------------------------------
def plot_map_2d(matfile):
  mat=mio.loadmat(matfile)
  XIsub=np.matrix(mat["XIsub"]/1000.0)
  YIsub=np.matrix(mat["YIsub"]/1000.0)
  ZIsub=np.matrix(mat["ZIsub"]/1000.0)

  return XIsub,YIsub,ZIsub
# -----------------------------------------------------------------------
def read_file(loc_filename):
  # Read location file
  if os.path.basename(loc_filename) == 'locations.dat' or os.path.basename(loc_filename) == 'relocations.dat':
    from locations_trigger import read_locs_from_file
    locs=read_locs_from_file(loc_filename)
  else:
    locs=[]
    with open(loc_filename,'r') as loc_file:
      loc_lines=loc_file.readlines()
      loc_file.close()
    for loc_line in loc_lines:
      word=loc_line.split()
      loc={}
      loc['o_time']=utcdatetime.UTCDateTime(word[0][:-1])
      loc['x_mean']=float(word[2])
      loc['y_mean']=float(word[5])
      loc['z_mean']=float(word[8])
      locs.append(loc)

  x,y,z,t,mag,amp=[],[],[],[],[],[]
  for loc in locs:
    x.append(loc['x_mean'])
    y.append(loc['y_mean'])
    z.append(-loc['z_mean'])
    t.append(float(loc['o_time']))
    if loc.has_key('max_trig'):
      amp.append(loc['max_trig'])
    if loc.has_key('ml'):
      if not math.isnan(loc['ml']):
        mag.append(loc['ml'])
      else:
        mag.append(0)

  xm,ym,zm=np.mean(x),np.mean(y),np.mean(z)
  et_x,et_y,et_z=np.std(x),np.std(y),np.std(z)
  swarm_ext=[xm-et_x-0.1,xm+et_x+0.1,ym-et_y-0.1,ym+et_y+0.1,zm-et_z-0.1,zm+et_z+0.1]
  print "Swarm: ",swarm_ext
  print "Nb evenements:", len(x)

  out_file=open('%s/out_swarm.dat'%os.path.dirname(loc_filename),'w')
  for loc in locs:
    if not xm-et_x-0.1 < loc['x_mean'] < xm+et_x+0.1 or not ym-et_y-0.1 < loc['y_mean'] < ym+et_y+0.1 or not zm-et_z-0.1 < -loc['z_mean'] < zm+et_z+0.1:
      out_file.write('%s, x= %s km, y= %s km, z= %s km\n'%(loc['o_time'].isoformat(), loc['x_mean'], loc['y_mean'], loc['z_mean']))

  out_file.close()

  return x,y,z,t,np.array(mag),swarm_ext
# -----------------------------------------------------------------------
def plot_seis_xy(x,y,ZIsub,ext,t,mag,sep=False):
  if mag.any():
    m=np.mean(mag)
    stdm=np.std(mag)
    plt.scatter(x,y,marker='o',c=mag,cmap=plt.cm.hot_r,s=10,vmin=m-stdm,vmax=m+stdm)
  elif sep:
    plt.plot(x,y,'bo',markersize=3)
    plt.plot(x[190:],y[190:],'yo',markersize=3)
  else:
    plt.scatter(x,y,c=t,cmap=plt.cm.gray)#,markersize=3)
  plt.imshow(ZIsub,cmap=plt.cm.gray,extent=(ext[0],ext[1],ext[2],ext[3]))
  plt.axis([ext[0],ext[1],ext[2],ext[3]])
  #plt.xlabel('x utm (km)',fontsize=18)
  plt.ylabel('y utm (km)',fontsize=18)
  #plt.colorbar()
# -----------------------------------------------------------------------
def plot_seis_xz(x,z,YIsub,ext,t,mag,sep=False):
  if type(YIsub) != np.ndarray:
    YIsub=np.array(YIsub)
  C=np.arange(ext[0],ext[1],.025)
  plt.plot(C,YIsub[:len(C)],'k-',alpha=.3)
  if mag.any():
    m=np.mean(mag)
    stdm=np.std(mag)
    plt.scatter(x,z,marker='o',c=mag,cmap=plt.cm.hot_r,s=10,vmin=m-stdm,vmax=m+stdm)
  elif sep:
    plt.plot(x,z,'bo',markersize=3)
    plt.plot(x[190:],z[190:],'yo',markersize=3)
  else:
    plt.scatter(x,z,c=t,cmap=plt.cm.gray)#,markersize=3)
  plt.axis([ext[0],ext[1],ext[4],ext[5]])
  plt.plot(range(YIsub.shape[0]),np.zeros(YIsub.shape[0]),'y--')
  plt.ylabel('z (km up)',fontsize=18)
  #plt.colorbar()
# -----------------------------------------------------------------------
def plot_seis_yz(y,z,XIsub,ext,t,mag,sep=False):
  if type(XIsub) != np.ndarray:
    XIsub=np.array(XIsub)
  C=np.arange(ext[2],ext[3],.025)
  plt.plot(C,XIsub[:len(C)],'k-',alpha=.3)
  if mag.any():
    m=np.mean(mag)
    stdm=np.std(mag)
    plt.scatter(y,z,marker='o',c=mag,cmap=plt.cm.hot_r,s=10,vmin=m-stdm,vmax=m+stdm)
  elif sep:
    plt.plot(y,z,'bo',markersize=3)
    plt.plot(y[190:],z[190:],'yo',markersize=3)
  else:
    plt.scatter(y,z,c=t,marker='o',cmap=plt.cm.gray)#,markersize=3)
  plt.axis([ext[2],ext[3],ext[4],ext[5]])
  plt.plot(C,np.zeros(len(C)),'y--')
  plt.xlabel('y utm (km)',fontsize=18)
  plt.ylabel('z (km up)',fontsize=18)
  #plt.colorbar()
# -----------------------------------------------------------------------
def slice_y(py,ext,YIsub,ZIsub):
  a=np.where(YIsub == py)[0][0,0]
  Z=np.array(ZIsub[a,:])
  return Z[0,:] 
# -----------------------------------------------------------------------
def slice_x(px,ext,XIsub,ZIsub):
  a=np.where(XIsub == px)[1][0,0]
  Z=np.array(ZIsub[:,a])
  Z=Z.T
  Z=Z[0,:]
  return Z[::-1]
# -----------------------------------------------------------------------
###### PLOT 3D - MAYAVI
def plot_3d_mayavi(matfile):
  XIsub,YIsub,ZIsub,ext=plot_map_3d(matfile)

  mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0),size=(1000,900))
  mlab.surf(YIsub,XIsub,ZIsub.transpose(),extent=ext,color=(0.8,0.8,0.8),opacity=.4)
  mlab.points3d(xsta,ysta,zsta,color=(1,0,0),scale_factor=0.1,mode='cube')

  #loc_filename="../data/14_oct_locpicks/loc_hyp/recap"
  loc_filename="../out/2010-10-14/loc/common.dat"
  x,y,z,t,swarm_ext=read_file(loc_filename)
  mlab.points3d(x,y,z,color=(0,0,1),scale_factor=0.1,opacity=1.)

  loc_filename="../out/2010-10-14/loc/new_common.dat"
  x,y,z,t,swarm_ext=read_file(loc_filename)
  mlab.points3d(x,y,z,color=(1,1,0),scale_factor=0.1,opacity=1.)

  mlab.axes(extent=ext,color=(0,0,0))
  mlab.outline(extent=ext,color=(0,0,0))


#plot_3d_mayavi(matfile)

###### PLOT 3D - MATPLOTLIB
#from mpl_toolkits.mplot3d import axes3d
#fig = plt.figure()
#fig.set_facecolor('white')
#ax = fig.add_subplot(111, projection='3d')
#print YIsub.shape,XIsub.shape,ZIsub.shape
#ax.plot_surface(XIsub,YIsub,ZIsub.transpose())


###### PLOT 2D - MATPLOTLIB
def plot_2d(loc_filename,xsta,ysta,zsta,size=False,sep=False,colsta=None):
  XIsub,YIsub,ZIsub=plot_map_2d(matfile)
  ZIsub=ZIsub[::-1]
  ext=[np.min(XIsub),np.max(XIsub),np.min(YIsub),np.max(YIsub),np.min(ZIsub),np.max(ZIsub)]

  x,y,z,t,mag,swarm_ext=read_file(loc_filename)
  if mag.any() and size:
    ev1=np.array([[x[i],y[i],z[i]] for i in np.where(mag >= 0)[0]])
    ev2=np.array([[x[i],y[i],z[i]] for i in np.where(mag >= 1)[0]])
    ev3=np.array([[x[i],y[i],z[i]] for i in np.where(mag >= 2)[0]])
  new_ext=[np.min(x),np.max(x),np.min(y),np.max(y),np.min(z),np.max(z)]

  for i in range(len(ext)):
    if i%2 == 0 and new_ext[i]-ext[i] < 0:
      ext[i]=new_ext[i]
    if i%2 !=0 and new_ext[i]-ext[i] > 0:
      ext[i]=new_ext[i]

  ext[4]=-1.0
  ext[5]=3.0

  fig=plt.figure(figsize=(9,10.5))
  fig.set_facecolor('white')
  import matplotlib.gridspec as gridspec
  G = gridspec.GridSpec(5,3)
  ax=fig.add_subplot(G[:3,:])
  plot_seis_xy(x,y,ZIsub,ext,t,mag,sep=sep)
  if colsta != 'r':
    ax.plot(xsta,ysta,marker='v',ls='',color=(1,0.8,.1),markersize=8)
  else:
    ax.plot(xsta,ysta,'rv',markersize=8)
  if mag.any() and size:
    ax.plot(x,y,'ko',markersize=3,markeredgecolor='k')
    if ev1.any() and ev1.all():
      ax.plot(ev1[:,0],ev1[:,1],'ro',markersize=3,markeredgecolor='r') 
    if ev2.any() and ev2.all():
      ax.plot(ev2[:,0],ev2[:,1],'yo',markersize=3,markeredgecolor='y')
    if ev3.any() and ev3.all():
      ax.plot(ev3[:,0],ev3[:,1],'wo',markersize=3,markeredgecolor='w')
  #if len(x) > 400:
    #ax.plot(x[410:],y[410:],'yo',markersize=3)
  ax.set_xticklabels('',minor=False)
  ax.tick_params(labelsize=18)
  ax.text(-0.15,0.98,"(i)",fontsize=16,transform=ax.transAxes)

  # add colorbar
  from mpl_toolkits.axes_grid1.inset_locator import inset_axes
  axins=inset_axes(ax,width="50%",height="5%",loc=2)
  plt.colorbar(cax=axins,orientation='horizontal')
  axins.xaxis.set_ticks_position("bottom")

  # plot magnitude histogram
  if mag.any(): 
    ax=fig.add_subplot(G[:1,2])
    a,b,c = ax.hist(mag,10,color='k',lw=1.5,alpha=.2)
    ax.set_xlabel('Magnitude',color='w')
    ax.set_xticks([-2,-1,0,1,2])
    ax.set_xticklabels([-2,-1,0,1,2],color='w')
    #pas=np.max(a)/100*20
    maxi=np.max(a)/10*10
    #yt=range(0,maxi+pas,pas)
    yt=range(0,maxi,50)
    ax.set_yticks(yt)
    ax.set_yticklabels(yt,color='w')
 
  py=7650
  Z=slice_y(py,ext,YIsub,ZIsub)
  ax=fig.add_subplot(G[3,:])
  plot_seis_xz(x,z,Z,ext,t,mag,sep=sep)
  if colsta != 'r':
    ax.plot(xsta,zsta,marker='v',ls='',color=(1,0.8,.1),markersize=8)
  else:
    ax.plot(xsta,zsta,'rv',markersize=8)
  if mag.any() and size:
    ax.plot(x,z,'ko',markersize=3,markeredgecolor='k')
    if ev1.any() and ev1.all():
      ax.plot(ev1[:,0],ev1[:,2],'ro',markersize=3,markeredgecolor='r')
    if ev2.any() and ev2.all():
      ax.plot(ev2[:,0],ev2[:,2],'yo',markersize=3,markeredgecolor='y')
    if ev3.any() and ev3.all():
      ax.plot(ev3[:,0],ev3[:,2],'wo',markersize=3,markeredgecolor='w')
  ax.set_yticks([-1,0,1,2,3])
  ax.text(-0.15,0.98,"(ii)",fontsize=16,transform=ax.transAxes)
  plt.text(ext[1]-6,ext[5]-0.7,'y=%d km'%py,fontsize=15)
  plt.text(ext[0]+0.5,ext[5]-0.7,'W',fontsize=15)
  plt.text(ext[1]-1,ext[5]-0.7,'E',fontsize=15)
  ax.set_xlabel('x utm (km)',fontsize=18,va='top')
  ax.xaxis.set_ticks_position('top')
  ax.tick_params(labelsize=18)

  px=366
  Z=slice_x(px,ext,XIsub,ZIsub)
  ax=fig.add_subplot(G[4,:])
  plot_seis_yz(y,z,Z,ext,t,mag,sep=sep)
  if colsta != 'r':
    ax.plot(ysta,zsta,marker='v',ls='',color=(1,0.8,.1),markersize=8)
  else:
    ax.plot(ysta,zsta,'rv',markersize=8)
  if mag.any() and size:
    ax.plot(y,z,'ko',markersize=3,markeredgecolor='k')
    if ev1.any() and ev1.all():
      ax.plot(ev1[:,1],ev1[:,2],'ro',markersize=3,markeredgecolor='r')
    if ev2.any() and ev2.all():
      ax.plot(ev2[:,1],ev2[:,2],'yo',markersize=3,markeredgecolor='y')
    if ev3.any() and ev3.all():
      ax.plot(ev3[:,1],ev3[:,2],'wo',markersize=3,markeredgecolor='w')
  ax.set_yticks([-1,0,1,2,3])
  plt.text(ext[3]-4,ext[5]-0.7,'x=%d km'%px,fontsize=15)
  plt.text(ext[2]+0.5,ext[5]-0.7,'S',fontsize=15)
  plt.text(ext[3]-0.5,ext[5]-0.7,'N',fontsize=15)
  ax.text(-0.15,0.98,"(iii)",fontsize=16,transform=ax.transAxes)
  ax.tick_params(labelsize=18)
  #fig.tight_layout()

  #fig.text(0.01,.95,"(a)",fontsize=18)

# -----------------------------------------------------------------------
if __name__ == '__main__' :

  # MAIN PROGRAM

  # Read station file
  sta_filename="../lib/coord_stations_piton"
  with open(sta_filename,'r') as sta:
    sta_lines=sta.readlines()
    sta.close()
  xsta,ysta,zsta,name=[],[],[],[]
  for sta_line in sta_lines:
    parts=sta_line.split()
    xsta.append(np.float(parts[3]))
    ysta.append(np.float(parts[4]))
    zsta.append(np.float(parts[6]))
    name.append(parts[1])

  matfile="../lib/MNT_PdF.mat"

  #crises=['2009-10-14','2009-10-18','2009-10-29','2009-11-05','2009-12-14','2009-12-29','2010-01-02','2010-09-23','2010-10-14','2010-12-09','2011-02-02']
  crises = ['2010-10-14']

  verbose=False
  format = 'png'

  for date in crises:

    print date
    if date == '2010-12-09':
      sep=True
    else:
      sep=False

    if date == '2010-10-14':
      loc_filename="../data/14_oct_locpicks/loc_hyp/recap"
      plot_2d(loc_filename,xsta,ysta,zsta,colsta='r')

      loc_filename="../out/%s/loc/locations.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,size=True)

      loc_filename="../out/%s/loc/relocations.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,size=True)

      loc_filename="../out/%s/loc/orig_reloc.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,colsta='r')

      loc_filename="../out/%s/loc/new_common.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,colsta='r')

      loc_filename="../out/%s/loc/reloc_common.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,colsta='r')

      loc_filename="../out/%s/loc/supp.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,colsta='r')

      loc_filename="../out/%s/loc/common.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,colsta='r')

    else:
      loc_filename="../out/%s/loc/locations.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,size=True,sep=sep)

      loc_filename="../out/%s/loc/relocations.dat"%date
      plot_2d(loc_filename,xsta,ysta,zsta,sep=sep,size=True)

    if verbose:
      plt.show()
    else:
      plt.close()
