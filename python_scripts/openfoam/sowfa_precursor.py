# Function for processing SOWFA source term data - in python.
#
# Origional matlab version by Matt Churchfield
# Converted to Python by Jen Annoni  - 2 Jan 2018
# Class structure implemented by Michael Lawson - Nov 2018

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
# from glob import glob
import read

class OpenFOAMInput(object):
    '''OpenFOAM input data object
    '''
    def __init__(self,input_file):
        self.data = {}
        self.file = os.path.abspath(input_file)
        self._read()


    def _read(self):
        '''
        Note that this function only reads scalar and boolean/string inputs from
        an OpenFOAM intput input file
        '''

        with open(self.file,'r') as fid:
            raw = fid.readlines()

        count = 0
        bloc_comment_test = False
        for i, line in enumerate(raw):

            if raw[i][0:2]=='/*':
                bloc_comment_test = True

            if bloc_comment_test is False:

                if raw[i].strip()[0:2]=='//' or raw[i].strip()[0:1]=='#': # Check if the string is a comment and skip line
                    pass

                elif len(raw[i].strip())==0: # Check if the string is empty and skip line
                    pass

                else:
                    tmp = raw[i].strip().rstrip().split()
                    try:
                        self.data[tmp[0].replace('"', '')] = np.float(tmp[1][:-1])
                    except:
                        self.data[tmp[0].replace('"', '')] = tmp[1][:-1]



            if raw[i][0:2]=='\*':
                bloc_comment_test = False

        print('Read SOWFA "setUp" file:\n\t' + self.file )

class Sim(object):
    def __init__(self, dir='./', log=None, time_dir=None, avg_time=None, avg_width=None, z_level=None):
        self.input = dict()
        self.output = dict()

        self.input['dir'] = dir
        self.input['avg_dir'] = os.path.join(dir,'postProcessing/averaging')
        self.input['time_dir'] = os.path.join(self.input['avg_dir'],time_dir)
        # self.input['times'] =  glob(os.path.join(self.input['avg_dir'],"*"))

        self.input['hLevelCellFile'] = os.path.join(self.input['time_dir'],'hLevelsCell')
        try:
            self.input['zCell'] = np.array(open(self.input['hLevelCellFile'],'r').read().split()).astype(np.float)
            print('Read "zCell" from:\n\t' + self.input['hLevelCellFile'])
        except:
            raise Exception('The file ' + self.input['hLevelCellFile'] + ' does not exist or is formatted incorrectly' )

        # Read setUp file
        setUpFile = os.path.join(dir,'setUp')
        self.setUp = OpenFOAMInput(input_file=setUpFile)
        self.input.update(self.setUp.data)

        # Get uStarMean from the log file
        self.input['logFile'] = os.path.join(dir,log)
        self._get_uStar()


        self.input['avg_time'] = avg_time
        self.input['avg_width'] = avg_width
        self.input['z_level'] = z_level

        # Check to make sure directories and files exist
        if os.path.isdir(self.input['avg_dir']) is True:
            pass
        else:
            raise Exception(self.input['avg_dir'] + ' does not exist')

        if os.path.isdir(self.input['time_dir']) is True:
            pass
        else:
            raise Exception(self.input['time_dir'] + ' does not exist')

    # def _create(self,dir=None,time_dir=None):


    def _get_uStar(self):
        with open(self.input['logFile']) as fid:
            raw = fid.readlines()
        for i in np.arange(-1,-100,-1):
            if np.shape(raw[i].split())[0] == 0:
                pass
            elif raw[i].split()[0] == 'uStarMean':
                self.input['uStarMean'] = np.float(raw[i].split()[2])
                break
        print('Set uStarMean = ' + str(self.input['uStarMean']) + ' from SOWFA ".log" file:\n\t' + self.input['logFile'])

    def theta_w_avg_cell(self):
        '''
        '''

        zCell = self.input['zCell']

        nzCell = len(zCell)

        # ni = len(self.input['times'])

        strData = os.path.join(self.input['time_dir'], 'Tw_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')

        t = np.array(data[0])
        dt = np.array(data[1])
        TwMean = data[list(range(2,len(data.columns)))]

        strData = os.path.join(self.input['time_dir'], 'q3_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')

        q3Mean = data[list(range(2,len(data.columns)))]

        # sum resolved and sub-grid stresses
        Twq3Mean = TwMean + q3Mean

        # perform average
        nt = len(t)
        nPts = len(Twq3Mean)

        # define start and end of averaging window
        tstart = np.max([self.input['avg_time'] - 0.5*self.input['avg_width'],t[0]])
        tend = np.min([tstart + self.input['avg_width'],t[len(t)-1]])

        # figure out the indices of the averaging window

        tidx1 = np.where(abs(t-tstart)==np.min(abs(t-tstart)))[0][0]
        tidx2 = np.where(abs(t-tend)==np.min(abs(t-tend)))[0][0]

        print('Time from ', t[tidx1], 'to ', t[tidx2])

        TwMeanAvg = np.zeros(len(Twq3Mean.columns))
        q3MeanAvg = np.zeros(len(Twq3Mean.columns))
        Twq3MeanAvg = np.zeros(len(Twq3Mean.columns))
        dtSum = 0
        for i in range(tidx1,tidx2):
            TwMeanAvg = TwMeanAvg + dt[i]*TwMean.iloc[i,:]
            q3MeanAvg = q3MeanAvg + dt[i]*q3Mean.iloc[i,:]
            Twq3MeanAvg = Twq3MeanAvg + dt[i]*Twq3Mean.iloc[i,:]
            dtSum = dtSum + dt[i]

        TwMeanAvg = TwMeanAvg/dtSum
        q3MeanAvg = q3MeanAvg/dtSum
        Twq3MeanAvg = Twq3MeanAvg/dtSum

        # find zi
        min_idx = np.where(Twq3MeanAvg==np.min(Twq3MeanAvg))
        ziAvg = zCell[min_idx[0][0]]

        # print('zi = ', ziAvg)


        plt.figure(num='theta_w_avg_cell')
        plt.plot(TwMeanAvg,zCell,'k--',linewidth=3,label=r'$\langle\theta w\rangle^{r}$')
        plt.plot(q3MeanAvg,zCell,'k:',linewidth=3,label=r'$\langle\theta w\rangle^{SGS}$')
        plt.plot(Twq3MeanAvg,zCell,'k',linewidth=3,label=r'$\langle\theta w\rangle^{Tot}$')
        plt.legend()
        plt.xlabel(r'$\langle \theta w \rangle$')
        plt.ylabel(r'$z/z_{i}$')

        plt.show()
        self.output['ziAvg'] = ziAvg

    def Umean_avg_nonnormalized(self):

        # Plots <U>

        zCell = self.input['zCell']

        nzCell = len(zCell)

        # U mean
        strData = os.path.join(self.input['time_dir'],  'U_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')

        t = np.array(data[0])
        dt = np.array(data[1])
        UMean = data[list(range(2,len(data.columns)))]

        # V mean
        strData = os.path.join(self.input['time_dir'], 'V_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')

        VMean = data[list(range(2,len(data.columns)))]

        # W mean
        strData = os.path.join(self.input['time_dir'], 'W_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')

        WMean = data[list(range(2,len(data.columns)))]

        nt = len(t)

        UVMean = np.sqrt((UMean**2 + VMean**2))

        # perform the average

        # define start and end of averaging window
        tstart = np.max([self.input['avg_time'] - 0.5*self.input['avg_width'],t[0]])
        tend = np.min([tstart + self.input['avg_width'],t[len(t)-1]])

        # figure out the indices of the averaging window
        tidx1 = np.where(abs(t-tstart)==np.min(abs(t-tstart)))[0][0]
        tidx2 = np.where(abs(t-tend)==np.min(abs(t-tend)))[0][0]

        UMeanAvg = np.zeros(len(UMean.columns))
        VMeanAvg = np.zeros(len(VMean.columns))
        WMeanAvg = np.zeros(len(WMean.columns))
        UVMeanAvg = np.zeros(len(UVMean.columns))

        dtSum = 0

        for i in range(tidx1,tidx2):
            UMeanAvg = UMeanAvg + dt[i]*UMean.iloc[i,:]
            VMeanAvg = VMeanAvg + dt[i]*VMean.iloc[i,:]
            WMeanAvg = WMeanAvg + dt[i]*WMean.iloc[i,:]
            UVMeanAvg = UVMeanAvg + dt[i]*UVMean.iloc[i,:]
            dtSum = dtSum + dt[i]

        UMeanAvg = UMeanAvg/dtSum
        VMeanAvg = VMeanAvg/dtSum
        WMeanAvg = WMeanAvg/dtSum
        UVMeanAvg = UVMeanAvg/dtSum

        # find d<U>dz at top of first cell for use later
        self.input['dUdz1Avg'] = (UMeanAvg.iloc[1]-UMeanAvg.iloc[0])/(zCell[1]-zCell[0])
        self.input['dvdz1Avg'] = (VMeanAvg.iloc[1]-UMeanAvg.iloc[0])/(zCell[1]-zCell[0])

        # phi_m
        z_f = 0.5*(zCell[1:] + zCell[0:-1])
        phi_m = (self.input['kappa']/self.input['uStarMean'])*z_f*((np.array(UVMeanAvg.iloc[1:])-np.array(UVMeanAvg.iloc[0:-1])))/(zCell[1:]-zCell[0:-1])


        # find wind direction
        windDir = np.array(np.arctan2(VMeanAvg,UMeanAvg))
        windDir = windDir*(180./np.pi)

        for i in range(len(windDir)):

            windDir[i] = 90.0 - windDir[i]

            if windDir[i] > 180.0:
                windDir[i] = windDir[i] - 180.0
            else:
                windDir[i] = windDir[i] + 180.0

            if (windDir[i] < 0.0):
                windDir[i] = windDir[i] + 360.0

        # find wind at various heights

        Uvec = np.zeros((len(self.input['heights']),3))
        Umag = np.zeros(len(self.input['heights']))

        self.input['dir'] = np.zeros(len(self.input['heights']))

        for i in range(len(self.input['heights'])):
            Uvec[i,0] = interp1d(zCell,UMeanAvg,kind='linear',fill_value='extrapolate')(self.input['heights'][i])
            Uvec[i,1] = interp1d(zCell,VMeanAvg,kind='linear',fill_value='extrapolate')(self.input['heights'][i])
            Uvec[i,2] = interp1d(zCell,WMeanAvg,kind='linear',fill_value='extrapolate')(self.input['heights'][i])

            Umag[i] = np.sqrt( Uvec[i,0]**2 + Uvec[i,1]**2 + Uvec[i,2]**2 )

            self.input['dir'][i] = interp1d(zCell,windDir,kind='linear',fill_value='extrapolate')(self.input['heights'][i])


        plt.figure(num='normalized u_mean_avg')
        plt.plot(phi_m,z_f/self.output['ziAvg'])
        plt.xlabel(r'$\phi_m$')
        plt.ylabel(r'$z/z_i$')
        plt.title('normalized u velocity')

        plt.figure(num='u_mean_avg, v_mean_avg, uv_mean_avg velocity')
        plt.plot(UMeanAvg,zCell,'k--',linewidth=3,label=r'$\langle U_x \rangle$')
        plt.plot(VMeanAvg,zCell,'k:',linewidth=3,label=r'$\langle U_y \rangle$')
        plt.plot(UVMeanAvg,zCell,'k-',linewidth=3,label=r'$\langle |U| \rangle$')
        plt.plot([0,16],[self.input['windHeight'],self.input['windHeight']],'r-')
        plt.legend()
        plt.xlabel(r'$\langle U_{i} \rangle$ (m/s)')
        plt.ylabel(r'$z$ (m)')
        plt.title('average velocities')

        plt.figure(num='w_mean velocity')
        plt.plot(WMeanAvg,zCell,'k-',linewidth=3)
        plt.xlabel(r'$\langle W \rangle$ (m/s)')
        plt.ylabel(r'$z$ (m)')
        plt.title('average w velocity')

        plt.figure(num='wind direction')
        plt.plot(windDir,zCell,'k-')
        plt.xlabel(r'wind direction ($^\circ$)')
        plt.ylabel(r'$z$ (m)')
        plt.title('wind direction')

        plt.figure(num='normalized velocities')
        plt.plot(UMeanAvg/self.input['U0Mag'],VMeanAvg/self.input['U0Mag'],'k-') #should this be wmag?
        plt.xlabel(r'$\langle U \rangle/U_{g}$')
        plt.ylabel(r'$\langle V \rangle/U_{g}$')
        plt.title('normalized velocities')

        plt.figure(num='Some turb stuff')
        plt.semilogx(zCell/self.input['z0'],self.input['kappa']*UVMeanAvg/self.input['uStarMean'],'k-')
        plt.xlabel(r'z/z_0')
        plt.ylabel(r'$(\kappa/u_{*}) \langle U \rangle$')
        plt.title(r'$\langle |U| \rangle$ law of the wall')

        self.input['Uvec'] = Uvec
        self.input['Umag'] = Umag

        # return self.input

    def Tmean_avg_nonnormalized(self):

        zCell = self.input['zCell']
        nzCell = len(zCell)

        strData = os.path.join(self.input['time_dir'], 'T_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')

        t = data[0]
        dt = data[1]
        TMean = data[list(range(2,len(data.columns)))]

        # get initial profile
        TMeanInitial = TMean.iloc[0,:]

        # perform average
        # define start and end of averaging window
        tstart = np.max([self.input['avg_time'] - 0.5*self.input['avg_width'],t[0]])
        tend = np.min([tstart + self.input['avg_width'],t[len(t)-1]])

        # figure out the indices of the averaging window
        tidx1 = np.where(abs(t-tstart)==np.min(abs(t-tstart)))[0][0]
        tidx2 = np.where(abs(t-tend)==np.min(abs(t-tend)))[0][0]

        TMeanAvg = np.zeros(len(TMean.columns))
        dtSum = 0

        for i in range(tidx1,tidx2):
            TMeanAvg = TMeanAvg + dt[i]*TMean.iloc[i,:]
            dtSum = dtSum + dt[i]

        TMeanAvg = TMeanAvg/dtSum

        # find wind at various heights

        Tmag = np.zeros(len(self.input['heights']))
        for i in range(len(self.input['heights'])):
            Tmag[i] = interp1d(zCell,TMeanAvg,kind='linear',fill_value='extrapolate')(self.input['heights'][i])
        self.input['Tmag'] = Tmag

        # plot results

        plt.figure(num='T vs zCell')
        plt.plot(TMeanAvg,zCell,'k-',label='Mean')
        plt.plot(TMeanInitial,zCell,'r:',label='Initial')
        plt.legend()
        plt.xlabel(r'$\langle \theta \rangle$ (K)')
        plt.ylabel(r'$z$ (m)')

        # return self.input

    def variances_avg_cell(self):

        # outputData = dict()

        # read in resolved stress data
        Umag = self.input['Umag']
        heights = self.input['heights']
        zCell = self.input['zCell']
        nzCell = len(zCell)

        # load uu
        strData = os.path.join(self.input['time_dir'], 'uu_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        t = data[0]
        dt = data[1]
        uuMean = data[list(range(2,len(data.columns)))]

        # load vv
        strData = os.path.join(self.input['time_dir'], 'vv_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        vvMean = data[list(range(2,len(data.columns)))]

        # load ww
        strData = os.path.join(self.input['time_dir'], 'ww_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        wwMean = data[list(range(2,len(data.columns)))]

        # load uv
        strData = os.path.join(self.input['time_dir'], 'uv_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        uvMean = data[list(range(2,len(data.columns)))]

        # load uw
        strData = os.path.join(self.input['time_dir'], 'uw_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        uwMean = data[list(range(2,len(data.columns)))]

        # load vw
        strData = os.path.join(self.input['time_dir'], 'vw_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        vwMean = data[list(range(2,len(data.columns)))]

        nt = len(t)

        # read in SFS stress data

        # load uu
        strData = os.path.join(self.input['time_dir'], 'R11_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        R11_Mean = data[list(range(2,len(data.columns)))]

        # load vv
        strData = os.path.join(self.input['time_dir'], 'R22_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        R22_Mean = data[list(range(2,len(data.columns)))]

        # load ww
        strData = os.path.join(self.input['time_dir'], 'R33_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        R33_Mean = data[list(range(2,len(data.columns)))]

        # load uv
        strData = os.path.join(self.input['time_dir'], 'R12_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        R12_Mean = data[list(range(2,len(data.columns)))]

        # load uw
        strData = os.path.join(self.input['time_dir'], 'R13_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        R13_Mean = data[list(range(2,len(data.columns)))]

        # load vw
        strData = os.path.join(self.input['time_dir'], 'R23_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        R23_Mean = data[list(range(2,len(data.columns)))]

        # perform average on resolved

        # perform average
        # define start and end of averaging window
        tstart = np.max([self.input['avg_time'] - 0.5*self.input['avg_width'],t[0]])
        tend = np.min([tstart + self.input['avg_width'],t[len(t)-1]])

        # figure out the indices of the averaging window
        tidx1 = np.where(abs(t-tstart)==np.min(abs(t-tstart)))[0][0]
        tidx2 = np.where(abs(t-tend)==np.min(abs(t-tend)))[0][0]

        uuMeanAvg = np.zeros(len(uuMean.columns))
        vvMeanAvg = np.zeros(len(vvMean.columns))
        wwMeanAvg = np.zeros(len(wwMean.columns))
        uvMeanAvg = np.zeros(len(uvMean.columns))
        uwMeanAvg = np.zeros(len(uwMean.columns))
        vwMeanAvg = np.zeros(len(vwMean.columns))

        R11_MeanAvg = np.zeros(len(R11_Mean.columns))
        R22_MeanAvg = np.zeros(len(R22_Mean.columns))
        R33_MeanAvg = np.zeros(len(R33_Mean.columns))
        R12_MeanAvg = np.zeros(len(R12_Mean.columns))
        R13_MeanAvg = np.zeros(len(R13_Mean.columns))
        R23_MeanAvg = np.zeros(len(R23_Mean.columns))

        dtSum = 0

        for i in range(tidx1,tidx2):
            uuMeanAvg = uuMeanAvg + dt[i]*uuMean.iloc[i,:]
            vvMeanAvg = vvMeanAvg + dt[i]*vvMean.iloc[i,:]
            wwMeanAvg = wwMeanAvg + dt[i]*wwMean.iloc[i,:]
            uvMeanAvg = uvMeanAvg + dt[i]*uvMean.iloc[i,:]
            uwMeanAvg = uwMeanAvg + dt[i]*uwMean.iloc[i,:]
            vwMeanAvg = vwMeanAvg + dt[i]*vwMean.iloc[i,:]

            R11_MeanAvg = R11_MeanAvg + dt[i]*R11_Mean.iloc[i,:]
            R22_MeanAvg = R22_MeanAvg + dt[i]*R22_Mean.iloc[i,:]
            R33_MeanAvg = R33_MeanAvg + dt[i]*R33_Mean.iloc[i,:]
            R12_MeanAvg = R12_MeanAvg + dt[i]*R12_Mean.iloc[i,:]
            R13_MeanAvg = R13_MeanAvg + dt[i]*R13_Mean.iloc[i,:]
            R23_MeanAvg = R23_MeanAvg + dt[i]*R23_Mean.iloc[i,:]

            dtSum = dtSum + dt[i]

        uuMeanAvg = uuMeanAvg/dtSum
        vvMeanAvg = vvMeanAvg/dtSum
        wwMeanAvg = wwMeanAvg/dtSum
        uvMeanAvg = uvMeanAvg/dtSum
        uwMeanAvg = uwMeanAvg/dtSum
        vwMeanAvg = vwMeanAvg/dtSum

        R11_MeanAvg = R11_MeanAvg/dtSum
        R22_MeanAvg = R22_MeanAvg/dtSum
        R33_MeanAvg = R33_MeanAvg/dtSum
        R12_MeanAvg = R12_MeanAvg/dtSum
        R13_MeanAvg = R13_MeanAvg/dtSum
        R23_MeanAvg = R23_MeanAvg/dtSum

        # sum the resolved and sfs stresses
        sum11_avg = uuMeanAvg + R11_MeanAvg
        sum22_avg = vvMeanAvg + R22_MeanAvg
        sum33_avg = wwMeanAvg + R33_MeanAvg
        sum12_avg = uvMeanAvg + R12_MeanAvg
        sum13_avg = uwMeanAvg + R13_MeanAvg
        sum23_avg = vwMeanAvg + R23_MeanAvg

        # find the turbulence intensity and tke at various heights
        windDir = np.arctan2(self.input['Uvec'][:,1],self.input['Uvec'][:,0])
        for i in range(len(windDir)):
            if windDir[i] < 0.0:
                windDir[i] = -1.0*windDir[i]

        TIxSum = np.zeros(len(heights))
        TIySum = np.zeros(len(heights))
        TIzSum = np.zeros(len(heights))

        tkeSum = np.zeros(len(heights))

        TIxyzSum = np.zeros(len(heights))

        TIdirSum = np.zeros(len(heights))

        TIxResolved = np.zeros(len(heights))
        TIyResolved = np.zeros(len(heights))
        TIzResolved = np.zeros(len(heights))

        tkeResolved = np.zeros(len(heights))

        TIxyzResolved = np.zeros(len(heights))

        TIdirResolved = np.zeros(len(heights))

        for i in range(len(heights)):

            TIxSum[i] = interp1d(zCell,np.sqrt(sum11_avg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])
            TIySum[i] = interp1d(zCell,np.sqrt(sum22_avg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])
            TIzSum[i] = interp1d(zCell,np.sqrt(sum33_avg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])

            tkeSum[i] = 0.5*interp1d(zCell,np.sqrt(sum11_avg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i]) + \
                            interp1d(zCell,np.sqrt(sum22_avg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i]) + \
                            interp1d(zCell,np.sqrt(sum33_avg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])

            TIxyzSum[i] = np.sqrt((2./3.)*tkeSum[i])/Umag[i]

            TIdirSum[i] = interp1d(zCell,sum11_avg,kind='linear',fill_value='extrapolate')(heights[i])*np.cos(windDir[i])*np.cos(windDir[i]) + \
                    2.0 * interp1d(zCell,sum12_avg,kind='linear',fill_value='extrapolate')(heights[i])*np.cos(windDir[i])*np.sin(windDir[i]) + \
                          interp1d(zCell,sum22_avg,kind='linear',fill_value='extrapolate')(heights[i])*np.sin(windDir[i])*np.sin(windDir[i])
            TIdirSum[i] = np.sqrt(TIdirSum[i])/Umag[i]

            TIxResolved[i] = interp1d(zCell,np.sqrt(uuMeanAvg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])
            TIyResolved[i] = interp1d(zCell,np.sqrt(vvMeanAvg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])
            TIzResolved[i] = interp1d(zCell,np.sqrt(wwMeanAvg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])

            tkeResolved[i] = 0.5*interp1d(zCell,np.sqrt(uuMeanAvg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i]) + \
                                 interp1d(zCell,np.sqrt(vvMeanAvg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i]) + \
                                 interp1d(zCell,np.sqrt(wwMeanAvg)/Umag[i],kind='linear',fill_value='extrapolate')(heights[i])

            TIxyzSum[i] = np.sqrt((2./3.)*tkeResolved[i])/Umag[i]

            TIdirResolved[i] = interp1d(zCell,uuMeanAvg,kind='linear',fill_value='extrapolate')(heights[i])*np.cos(windDir[i])*np.cos(windDir[i]) + \
                         2.0 * interp1d(zCell,vvMeanAvg,kind='linear',fill_value='extrapolate')(heights[i])*np.cos(windDir[i])*np.sin(windDir[i]) + \
                               interp1d(zCell,wwMeanAvg,kind='linear',fill_value='extrapolate')(heights[i])*np.sin(windDir[i])*np.sin(windDir[i])
            TIdirResolved[i] = np.sqrt(TIdirSum[i])/Umag[i]


        plt.figure(num='uu vs zCell')
        plt.plot(sum11_avg,zCell,'k--',label=r'$\langle uu \rangle$')
        plt.plot(sum22_avg,zCell,'k:',label=r'$\langle vv \rangle$')
        plt.plot(sum33_avg,zCell,'k-',label=r'$\langle ww \rangle$')
        plt.plot(uuMeanAvg,zCell,'b--')
        plt.plot(vvMeanAvg,zCell,'b:')
        plt.plot(wwMeanAvg,zCell,'b-')
        plt.plot(R11_MeanAvg,zCell,'r--')
        plt.plot(R22_MeanAvg,zCell,'r:')
        plt.plot(R33_MeanAvg,zCell,'r-')

        plt.xlabel(r'$\langle u_i u_j \rangle^{r} $ (m$^2$/s$^2$)')
        plt.ylabel(r'$z/z_i$')

        plt.legend()


        plt.figure(num='uu vs zCell (2)')
        plt.plot(uuMeanAvg,zCell,'k--',label=r'$\langle uu \rangle$')
        plt.plot(vvMeanAvg,zCell,'k:',label=r'$\langle vv \rangle$')
        plt.plot(wwMeanAvg,zCell,'k-',label=r'$\langle ww \rangle$')

        plt.xlabel(r'$\langle u_{i} u_{j} \rangle^r$ (m$^2$/s$^2$)')
        plt.ylabel(r'$z/z_i$')

        plt.legend()

        self.output['TIxResolved'] = TIxResolved
        self.output['TIyResolved'] = TIyResolved
        self.output['TIzResolved'] = TIzResolved
        self.output['TIxyzResolved'] = TIxyzResolved
        self.output['TIdirResolved'] = TIdirResolved
        self.output['tkeResolved'] = tkeResolved
        self.output['TIxSum'] = TIxSum
        self.output['TIySum'] = TIySum
        self.output['TIzSum'] = TIzSum
        self.output['TIxyzSum'] = TIxyzSum
        self.output['TIdirSum'] = TIdirSum
        self.output['tkeSum'] = tkeSum


        # return self.input

    def Umean_h(self):

        # plots <U> and <V> at a specified height vs. time

        zCell = self.input['zCell']

        zindex = np.where(abs(zCell-self.input['z_level'])==np.min(abs(zCell-self.input['z_level'])))[0][0]

        if zCell[zindex] < self.input['z_level']:
            zindex1 = zindex
            zindex2 = zindex+1
        else:
            zindex1 = zindex-1
            zindex2 = zindex

        strData = os.path.join(self.input['time_dir'], 'U_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        t = data[0]
        dt = data[1]
        Umean = np.array(data[list(range(2,len(data.columns)))])

        strData = os.path.join(self.input['time_dir'], 'V_mean')
        data = pd.read_csv(strData,header=None,delimiter=' ')
        Vmean = np.array(data[list(range(2,len(data.columns)))])

        nt = len(t)

        # find <U> and <V> at a specified height
        Umeanz = Umean[:,zindex1] + (self.input['z_level']-zCell[zindex1])*(Umean[:,zindex2]-Umean[:,zindex1])/(zCell[zindex2]-zCell[zindex1])
        Vmeanz = Vmean[:,zindex1] + (self.input['z_level']-zCell[zindex1])*(Vmean[:,zindex2]-Vmean[:,zindex1])/(zCell[zindex2]-zCell[zindex1])

        # plot stuff
        plt.figure(num='t vs Umeanz/U0Mag')
        plt.plot(t,Umeanz/self.input['U0Mag'],'k')
        plt.xlabel('t (s)')
        plt.ylabel(r'$\langle U \rangle /U_g$')


        plt.figure(num='t vs Vmeanz/U0Mag')
        plt.plot(t,Vmeanz/self.input['U0Mag'],'k')
        plt.xlabel('t (s)')
        plt.ylabel(r'$\langle V \rangle /U_g$')
