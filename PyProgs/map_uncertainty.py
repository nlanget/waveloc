import sys,h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def mapping():

  from NllGridLib import read_stations_file
  stations = read_stations_file('/home/nadege/waveloc/lib/coord_stations_ijen_utm')

  nx,ny,nz = 41,41,11
  dx,dy,dz = 0.5,0.5,0.5
  x_orig,y_orig,z_orig = -10,-10,-3
  x = np.arange(x_orig,x_orig+dx*nx,dx)
  y = np.arange(y_orig,y_orig+dy*ny,dy)
  z = np.arange(z_orig,z_orig+dz*nz,dz)

  xsta = np.array([np.argmin(np.abs(stations[key]['x']-x)) for key in sorted(stations)])
  ysta = np.array([np.argmin(np.abs(stations[key]['y']-y)) for key in sorted(stations)])
  zsta = np.array([np.argmin(np.abs(-stations[key]['elev']-z)) for key in sorted(stations)])

  f = h5py.File('/home/nadege/waveloc/out/TEST_Dirac/incertitude.hdf5','r')
  #f = h5py.File('/home/nadege/waveloc/out/TEST_Dirac/uncertainty_ok.hdf5','r')
  unc_grid_x = f['unc_grid_x']
  unc_grid_y = f['unc_grid_y']
  unc_grid_z = f['unc_grid_z']

  zeros = np.where(unc_grid_x[:]==0)[0]
  zeros = np.append(zeros,np.where(unc_grid_y[:]==0)[0])
  zeros = np.append(zeros,np.where(unc_grid_z[:]==0)[0])
  zeros = np.unique(zeros)
  print zeros

  argsort = np.argsort(unc_grid_x[:])
  print argsort[-80:]
  print np.array(unc_grid_x)[argsort[-80:]]
  sys.exit()

  #f = h5py.File('/home/nadege/waveloc/out/TEST_Dirac/erreur.hdf5','r')
  #unc_grid_x = f['err_grid_x']
  #unc_grid_y = f['err_grid_y']
  #unc_grid_z = f['err_grid_z']
  grid_3D_x = unc_grid_x[:].reshape(nx,ny,nz)
  grid_3D_y = unc_grid_y[:].reshape(nx,ny,nz)
  grid_3D_z = unc_grid_z[:].reshape(nx,ny,nz)

  print "x",np.max(grid_3D_x),np.min(grid_3D_x),np.mean(grid_3D_x)
  print "y",np.max(grid_3D_y),np.min(grid_3D_y),np.mean(grid_3D_y)
  print "z",np.max(grid_3D_z),np.min(grid_3D_z),np.mean(grid_3D_z)
  vmin, vmax = 0, 2.5
  print "x",np.argmax(grid_3D_x),np.argmin(grid_3D_x)
  print "y",np.argmax(grid_3D_y),np.argmin(grid_3D_y)
  print "z",np.argmax(grid_3D_z),np.argmin(grid_3D_z)

  ix_true,iy_true,iz_true = 25,25,5

  x_xy_cut = grid_3D_x[:,:,iz_true]
  x_xz_cut = grid_3D_x[:,iy_true,:]
  x_yz_cut = grid_3D_x[ix_true,:,:]

  y_xy_cut = grid_3D_y[:,:,iz_true]
  y_xz_cut = grid_3D_y[:,iy_true,:]
  y_yz_cut = grid_3D_y[ix_true,:,:]

  z_xy_cut = grid_3D_z[:,:,iz_true]
  z_xz_cut = grid_3D_z[:,iy_true,:]
  z_yz_cut = grid_3D_z[ix_true,:,:]

  xlab = np.array(np.append([0],x[0::10]),dtype=int)
  ylab = np.array(np.append(y[0::10],[0]),dtype=int)
  #xlab = [0,-10,-7.5,-5,-2.5,0,2.5,5,0]
  #ylab = xlab
  zlab = [0,3,2,1,0,-1,-2]

  G = gridspec.GridSpec(3,3)
  fig = plt.figure(figsize=(10,9))
  fig.set_facecolor('white')
  ax = fig.add_subplot(G[0,0])
  ax.plot(xsta,np.abs(ysta-len(y)),'kv')
  cax = ax.imshow(x_xy_cut.T[::-1],vmin=0,vmax=vmax)
  ax.set_xlabel('x utm (km)')
  ax.set_ylabel('y utm (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(ylab[::-1])
  ax.set_title('Uncertainty on x')
  pos = list(ax.get_position().bounds)
  fig.text(pos[0]-0.08,pos[1]+pos[3], 'z=%.1f km'%-z[iz_true], fontsize=12)
  cbar = fig.colorbar(cax,ticks=[vmin,(vmin+vmax)/2.,vmax])
  cbar.ax.set_yticklabels([vmin,(vmin+vmax)/2.,'> %.1f'%vmax])

  ax = fig.add_subplot(G[1,0])
  ax.plot(xsta,zsta,'kv')
  ax.imshow(x_xz_cut.T,vmin=0,vmax=vmax)
  ax.set_xlabel('x utm (km)')
  ax.set_ylabel('z (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(zlab)
  ax.set_title('Uncertainty on x')
  pos = list(ax.get_position().bounds)
  fig.text(pos[0]-0.08,pos[1]+pos[3], 'y=%.1f km'%y[iy_true], fontsize=12)

  ax = fig.add_subplot(G[2,0])
  ax.plot(ysta,zsta,'kv')
  ax.imshow(x_yz_cut.T,vmin=0,vmax=vmax)
  ax.set_xlabel('y utm (km)')
  ax.set_ylabel('z (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(zlab)
  ax.set_title('Uncertainty on x')
  pos = list(ax.get_position().bounds)
  fig.text(pos[0]-0.08,pos[1]+pos[3], 'x=%.1f km'%x[ix_true], fontsize=12)

  ax = fig.add_subplot(G[0,1])
  ax.plot(xsta,np.abs(ysta-len(y)),'kv')
  cax = ax.imshow(y_xy_cut.T[::-1],vmin=0,vmax=vmax)
  ax.set_xlabel('x utm (km)')
  #ax.set_ylabel('y utm (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(ylab[::-1])
  ax.set_title('Uncertainty on y')
  cbar = fig.colorbar(cax,ticks=[vmin,(vmin+vmax)/2.,vmax])
  cbar.ax.set_yticklabels([vmin,(vmin+vmax)/2.,'> %.1f'%vmax])

  ax = fig.add_subplot(G[1,1])
  ax.plot(xsta,zsta,'kv')
  ax.imshow(y_xz_cut.T,vmin=0,vmax=vmax)
  ax.set_xlabel('x utm (km)')
  ax.set_ylabel('z (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(zlab)
  ax.set_title('Uncertainty on y')

  ax = fig.add_subplot(G[2,1])
  ax.plot(ysta,zsta,'kv')
  ax.imshow(y_yz_cut.T,vmin=0,vmax=vmax)
  ax.set_xlabel('y utm (km)')
  ax.set_ylabel('z (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(zlab)
  ax.set_title('Uncertainty on y')

  ax = fig.add_subplot(G[0,2])
  ax.plot(xsta,np.abs(ysta-len(y)),'kv')
  cax = ax.imshow(z_xy_cut.T[::-1],vmin=0,vmax=vmax)
  ax.set_xlabel('x utm (km)')
  #ax.set_ylabel('y utm (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(ylab[::-1])
  ax.set_title('Uncertainty on z')
  cbar = fig.colorbar(cax,ticks=[vmin,(vmin+vmax)/2.,vmax])
  cbar.ax.set_yticklabels([vmin,(vmin+vmax)/2.,'> %.1f'%vmax])

  ax = fig.add_subplot(G[1,2])
  ax.plot(xsta,zsta,'kv')
  ax.imshow(z_xz_cut.T,vmin=0,vmax=vmax)
  ax.set_xlabel('x utm (km)')
  ax.set_ylabel('z (km)') 
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(zlab)
  ax.set_title('Uncertainty on z')

  ax = fig.add_subplot(G[2,2])
  ax.plot(ysta,zsta,'kv')
  ax.imshow(z_yz_cut.T,vmin=0,vmax=vmax)
  ax.set_xlabel('y utm (km)')
  ax.set_ylabel('z (km)')
  ax.set_xticklabels(xlab)
  ax.set_yticklabels(zlab)
  ax.set_title('Uncertainty on z')

  plt.savefig('/home/nadege/Desktop/unc.png')


#  fig = plt.figure()
#  fig.set_facecolor('white')
#  ax = fig.add_subplot(311)
#  ax.plot(xsta,ysta,'kv')
#  ax.set_xlim([0,len(x)])
#  ax.set_ylim([0,len(y)])
#  ax = fig.add_subplot(312)
#  ax.plot(xsta,zsta,'kv')
#  ax.set_xlim([0,len(x)])
#  ax.set_ylim([len(z),0])
#  ax = fig.add_subplot(313)
#  ax.plot(ysta,zsta,'kv')
#  ax.set_xlim([0,len(y)])
#  ax.set_ylim([len(z),0])

  plt.show()



if __name__ == '__main__':
   mapping()