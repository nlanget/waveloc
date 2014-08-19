#!/usr/bin/env python
# encoding: utf-8

import glob,os
import numpy as np
import matplotlib.pyplot as plt
from NllGridLib import qd_read_hyp_file, qd_read_picks_from_hyp_file
from obspy.core import utcdatetime
import matplotlib.mlab as mlab
from scipy.interpolate import UnivariateSpline
from locations_trigger import read_locs_from_file, read_header_from_file

def histo(locdirs,loca):
  path="/home/nadege/waveloc"

  # Manual locations
  hyp_path="%s/data/14_oct_locpicks/loc_hyp"%path
  hyp_files=glob.glob("%s/piton*.hyp"%hyp_path)
  hyp_parameters=[]
  recap=open('/home/nadege/waveloc/data/14_oct_locpicks/loc_hyp/recap','w')
  for hyp_file in hyp_files:
    (otime,hypo_x,sigma_x,hypo_y,sigma_y,hypo_z,sigma_z)=qd_read_hyp_file(hyp_file)
    phase_dict=qd_read_picks_from_hyp_file(hyp_file)
    hyp_parameters.append(((otime,hypo_x,sigma_x,hypo_y,sigma_y,hypo_z,sigma_z),phase_dict))
    recap.write('%s, x= %s km, y= %s km, z= %s km\n'%(otime, hypo_x, hypo_y, hypo_z))
  number_hyp=len(hyp_files)
  recap.close()

  dx,dy,dz,dt,d,inc={},{},{},{},{},{}
  spec=[]
  a,b,mt,mx,my,mz,md,minc=[],[],[],[],[],[],[],[]
  number_of_locations=[]
  common=np.zeros(len(locdirs))
  i=0

  resdir = '/home/nadege/waveloc/out/Piton'
  for locdir in locdirs:
    com_file=open('%s/%s/common.dat'%(resdir,locdir),'w')
    supp=open('%s/%s/supp.dat'%(resdir,locdir),'w')
    dx[i],dy[i],dz[i],dt[i],d[i],inc[i]=[],[],[],[],[],[]
    locpath="%s/%s"%(resdir,locdir)
    loc_filename="%s/%s"%(locpath,loca)

    #if os.path.isfile("%s/relocations.dat"%locpath):
      #loc_filename="%s/relocations.dat"%locpath

    # Read location file
    locs=read_locs_from_file(loc_filename)
    head=read_header_from_file(loc_filename,{})
    number_of_locations.append(len(locs))

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
        common[i]+=1
        #com_file.write('%s, x= %s km, y= %s km, z= %s km\n'%(stack_time, stack_x, stack_y, stack_z))
        com_file.write('%s, x= %s km, y= %s km, z= %s km\n'%(otime, hypo_x, hypo_y, hypo_z))
        #com_file.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(stack_max,stack_time.isoformat(),stack_time_err_left, stack_time_err_right ,stack_x,stack_x_err,stack_y,stack_y_err,stack_z,stack_z_err))
        dt[i].append(stack_time-otime)
        dz[i].append(-stack_z+hypo_z)
        dy[i].append(stack_y-hypo_y)
        dx[i].append(stack_x-hypo_x)
        dauto=np.sqrt(stack_x**2+stack_y**2+stack_z**2)
        dman=np.sqrt(hypo_x**2+hypo_y**2+hypo_z**2)
        dist=np.sqrt((stack_x-hypo_x)**2+(stack_y-hypo_y)**2+(-stack_z+hypo_z)**2)
        d[i].append(dist)
        inc[i].append(1/np.sqrt(stack_x**2+stack_y**2+stack_z**2)*np.sqrt((stack_x**2*stack_x_err**2+stack_y**2*stack_y_err**2+stack_z**2*stack_y_err**2)))

      else:
        #supp.write("Max = %.2f, %s - %.2f s + %.2f s, x= %.4f pm %.4f km, y= %.4f pm %.4f km, z= %.4f pm %.4f km\n"%(stack_max,stack_time.isoformat(),stack_time_err_left, stack_time_err_right ,stack_x,stack_x_err,stack_y,stack_y_err,stack_z,stack_z_err))
        supp.write('%s, x= %s km, y= %s km, z= %s km\n'%(stack_time.isoformat(), stack_x, stack_y, stack_z))

    mt.append(np.mean(dt[i]))
    mx.append(np.mean(dx[i]))
    my.append(np.mean(dy[i]))
    mz.append(np.mean(dz[i]))
    md.append(np.mean(d[i]))
    minc.append(np.mean(inc[i]))

    comment=''
    if head['krec'] == 'True':
      comment='rec'
    if head['kderiv'] == 'True':
      comment='%s deriv'%comment
    if head['gauss'] == 'True':
      comment='%s gauss'%comment
    if head['reloc'] == 'True':
      comment='%s reloc'%comment

    spec.append(r'%.0f-%.0fHz %s ($K_t$=%s)'%(head['c1'],head['c2'],comment,head['loclevel']))

    i=i+1

  lab={}
  lab['time'],lab['x'],lab['y'],lab['z']=[],[],[],[]
  for i in range(len(d)):
    l='%s \n(W: %d - C: %d - M: %.2f)'%(spec[i],number_of_locations[i],common[i],mt[i])
    lab['time'].append(l)
    l='%s \n(W: %d - C: %d - M: %.2f)'%(spec[i],number_of_locations[i],common[i],mx[i])
    lab['x'].append(l)
    l='%s \n(W: %d - C: %d - M: %.2f)'%(spec[i],number_of_locations[i],common[i],my[i])
    lab['y'].append(l)
    l='%s \n(W: %d - C: %d - M: %.2f)'%(spec[i],number_of_locations[i],common[i],mz[i])
    lab['z'].append(l)

    i=i+1

  plot(number_of_locations,common,number_hyp,mt,dt,lab,'time')
  plot(number_of_locations,common,number_hyp,mx,dx,lab,'x')
  plot(number_of_locations,common,number_hyp,my,dy,lab,'y')
  plot(number_of_locations,common,number_hyp,mz,dz,lab,'z')

  com_file.close()
  supp.close()

# -------------------------------------------------------------------------- #
def spline(bins,n):
    col=['b','g','r','c','m','y']
    for i in range(len(n)):
      c=col[i]
      if i==0:
        s1=UnivariateSpline(bins[0:-1],n[i])
        s2=UnivariateSpline(bins,s1(bins))
        y=s2(bins)
      else:
        s1=UnivariateSpline(bins[0:-1],n[i])#,s=50)
        s2=UnivariateSpline(bins,s1(bins))
        s3=UnivariateSpline(bins,s2(bins))
        y=s3(bins)
      maxi=np.argmax(y)
      #for j in range(len(y)):
        #if j < maxi and y[j] > y[j+1]:
          #if j == 0:
            #y[j]=0
          #else:
            #y[j]=(y[j-1]+y[j+1])/2
        #elif j > maxi and y[j] > y[j-1]:
          #if j == len(y)-1:
            #y[j]=0
          #else:
            #y[j]=(y[j-1]+y[j+1])/2
      #plt.plot(bins[:-1],n[i],color=c)
      plt.plot(bins,y,color=c,linewidth=1.5)
      #plt.fill(bins,y,color=c,edgecolor='k',alpha=0.4)

# -------------------------------------------------------------------------- #
def plot(number_of_locations,common,number_hyp,m,d,lab,comp):

  fig = plt.figure(figsize=(11,9))
  fig.set_facecolor('white')

  Label=lab[comp]
  ax=fig.add_subplot(111)#,title='$\Delta$%s'%comp)
  if comp == 'time':
    let = '(a)'
  elif comp == 'x':
    let = '(b)'
  elif comp == 'y':
    let = '(c)'
  else:
    let = '(d)'

  if comp ==  'time':
    int=[-1,1]
    n, bins, patches = ax.hist(d.values(),50,range=None,align='left',histtype='stepfilled',alpha=.2,label=Label)
    #plt.figtext(0.1,0.6,"manual picks: %d"%number_hyp)
    unit='s' 

  else:
    int=[-3,3]
    n, bins, patches = ax.hist(d.values(),100,range=None,align='left',histtype='stepfilled',alpha=.2,label=Label)
    #plt.figtext(0.1,0.6,"manual picks: %d"%number_hyp)
    unit='km'

  ax.text(-0.12,0.99,let,fontsize=24,transform=ax.transAxes)
  ax.legend(loc=2,prop={'size':'x-large'})
  ax.axis(int+[0,np.max(n)+3])
  spline(bins,n)
  for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
  for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(20)

  ax.set_xlabel('$\Delta$%s (%s)'%(comp,unit),fontsize=24)
  ax.set_ylabel('Number of common events',fontsize=24)

  plt.savefig('/home/nadege/Desktop/histo_%s.png'%comp)

def plot_kurto(outdir):
  from correlation import BinaryFile
  kurto_file="%s/kurto"%outdir
  a=BinaryFile(kurto_file)
  param=a.read_binary_file()

  s1,s2=0,0
  
  for sta in sorted(param):
    fig=plt.figure()
    fig.set_facecolor('white')
    plt.hist(param[sta][:,0],50)
    plt.title(sta)

    fig=plt.figure()
    fig.set_facecolor('white')
    plt.hist(param[sta][:,1],50)
    plt.title(sta)
    #plt.show()

    print sta,np.mean(param[sta][:,0]),np.mean(param[sta][:,1])
    s1+=np.mean(param[sta][:,0])
    s2+=np.mean(param[sta][:,1])

  print s1/len(sorted(param)),s2/len(sorted(param))

#########################################################################
#locdirs = ["2010-10-14_20-35_grad/loc","2010-10-14_20-35_rec_1s_grad/loc","2010-10-14/loc"]
#locdirs = ["2010-10-14_4-10/loc","2010-10-14_20-35/loc","2010-10-14_20-35_grad/loc"]
#locdirs = ["2010-10-14_20-35_rec_3s/loc","2010-10-14_20-35_rec_1s/loc","2010-10-14_20-35_rec_05s/loc"]

locdirs = ["2010-10-14_4-10_rec/loc","2010-10-14_4-10/loc"]
histo(locdirs,'locations.dat')

#outdir='/home/nadege/waveloc/out/2010-10-14_4-10_rec'
#plot_kurto(outdir)

#plt.show()

