import os,glob,sys
import pandas as pd
import numpy as np
from obspy.core import read,utcdatetime
from correlation import correlate,BinaryFile

def do_ijen_correlation():
  filename = '/home/nadege/Desktop/IJEN_catalogues/test_extraction.csv'
  datapath = '/home/nadege/Desktop/NEW_CLASS/Cat_POS/POS'
 
  tm, tp = 5,20  # time before and after the origin time
  df = pd.read_csv(filename)
  coeff,delay = [],[]
  for idate,date in enumerate(df.Date):
    date = "%s_%s"%(str(date)[:8],str(date)[8:])
    date_utc = utcdatetime.UTCDateTime(df.Date_UTC[idate])
    files = glob.glob(os.path.join(datapath,'*Z*%s*'%date))
    if len(files) > 0:
      file = files[0]
      st1 = read(file,starttime=date_utc-tm,endtime=date_utc+tp)
      st1.filter('bandpass',freqmin=1,freqmax=10)
      tr = st1[0]
      x = tr.data
      dt = tr.stats.delta

      sub_coeff,sub_delay = [],[]
      for id,d in enumerate(df.Date[idate+1:]):
        d = "%s_%s"%(str(d)[:8],str(d)[8:])
        d_utc = utcdatetime.UTCDateTime(df.Date_UTC[idate+1+id])
        fs = glob.glob(os.path.join(datapath,'*Z*%s*'%d))
        if len(fs) > 0:
          f = fs[0]
          st2 = read(f,starttime=d_utc-tm,endtime=d_utc+tp)
          st2.filter('bandpass',freqmin=1,freqmax=10)
          tr = st2[0]
          y = tr.data

          tau,value = correlate(x,y,dt,False,'t')
          sub_coeff.append(value)
          sub_delay.append(tau)

          #if value > 0.75:
          #  print file,f
          #  print tau,value
          #  st1.plot()
          #  st2.plot()

      coeff.append(sub_coeff)
      delay.append(sub_delay)

  a = BinaryFile('%s/MyClass_Coeff'%os.path.dirname(filename))
  a.write_binary_file(coeff)
  b = BinaryFile('%s/MyClass_Tau'%os.path.dirname(filename))
  b.write_binary_file(delay)


def interp_corr():
  import matplotlib.pyplot as plt
  dir = '/home/nadege/Desktop/IJEN_catalogues'
  a = BinaryFile('%s/MyClass_Coeff'%dir)
  coeff = a.read_binary_file()
  b = BinaryFile('%s/MyClass_Tau'%dir)
  delay = b.read_binary_file()

  dim = len(coeff[0])+1
  for i in range(dim):
    if len(coeff[i]) != dim:
      z = list(np.zeros(dim-len(coeff[i])))
      coeff[i] = z + coeff[i]
      delay[i] = z + delay[i]
  coeff = np.array(coeff)
  delay = np.array(delay)

  threshold = .75
  array = np.where(coeff > threshold)

  df = pd.read_csv('%s/test_extraction.csv'%dir)
  datapath = '/home/nadege/Desktop/NEW_CLASS/Cat_POS/POS'
  for i in range(array[0].shape[0]):
    ev1 = array[0][i]
    d1 = df.Date[ev1]
    d1 = "%s_%s"%(str(d1)[:8],str(d1)[8:])

    ev2 = array[1][i]
    d2 = df.Date[ev2]
    d2 = "%s_%s"%(str(d2)[:8],str(d2)[8:])
    
    f1 = glob.glob(os.path.join(datapath,'*Z*%s*'%d1))[0]
    f2 = glob.glob(os.path.join(datapath,'*Z*%s*'%d2))[0]

    tdeb = utcdatetime.UTCDateTime(df.Date_UTC[ev1])-5
    tfin = utcdatetime.UTCDateTime(df.Date_UTC[ev1])+20
    st1 = read(f1,starttime=tdeb,endtime=tfin)
    st1.filter('bandpass',freqmin=1,freqmax=10)
    tr = st1[0]
    data_1 = tr.data
    dt = tr.stats.delta

    tdeb = utcdatetime.UTCDateTime(df.Date_UTC[ev2])-5
    tfin = utcdatetime.UTCDateTime(df.Date_UTC[ev2])+20
    st2 = read(f2,starttime=tdeb,endtime=tfin)
    st2.filter('bandpass',freqmin=1,freqmax=10)
    tr = st2[0]
    data_2 = tr.data


    print df.Date_UTC[ev1],df.Date_UTC[ev2]
    print "Correlation coefficient :",coeff[ev1,ev2]
    print "Delay :",delay[ev1,ev2]
    print "Types :",df.Type[ev1],df.Type[ev2]
    fig = plt.figure()
    fig.set_facecolor('white')
    ax1 = fig.add_subplot(311)
    ax1.plot(data_1,'k')
    ax2 = fig.add_subplot(312)
    ax2.plot(data_2,'k')
    ax3 = fig.add_subplot(313)
    ax3.plot(data_1/np.max(data_1),'k')
    num = int(1./dt*delay[ev1,ev2])
    ax3.plot(range(num,len(data_1)+num),data_2/np.max(data_2),'r')
    plt.show()


if __name__ == '__main__':
  #do_ijen_correlation()
  interp_corr()
