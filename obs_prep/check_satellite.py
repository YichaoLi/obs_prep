from astropy import units as u
from skyfield.api import Topos, load
import numpy as np
from astropy.time import Time, TimezoneInfo
from astropy.coordinates import AltAz
import copy

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.dates as mdates

#from Queue import Queue
#import threading
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial

import time

_c_list = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
          "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]

_gps_tles = [ 'gps-ops', 'glo-ops', 'beidou', 'galileo', 'sbas', ] #'active', 
_comms_tles = ['geo',    'iridium', 'iridium-NEXT', 'nnss']



def get_sat_tles(sattype='geo', reload=False, source_url=None):
    '''
    Collect the satellite TLEs that fall in the L-Band.
    Parameters
    ----------
    geostationary: boolean
        Include geostationary staellites or not.
    Returns
    -------
    sats: dictionary
        Satellite TLEs ready for consumption by Skyfield.
    '''

    if source_url is None:
        source_url = 'https://www.celestrak.com/NORAD/elements/'
    #active = get_active_tles(source_url=source_url)
    print('load sat catalogue from %s %s.txt'%(source_url, sattype))
    LbandSats = [sattype + '.txt', ]
    sats = {}
    for i in range(len(LbandSats)):
        _sats = load.tle(source_url+LbandSats[i], reload=reload)
        for xx in list(_sats.keys()):
            if not isinstance(xx, int): # and xx in active.keys():
                sats[xx] = _sats[xx]

    return sats

def get_active_tles(source_url=None):

    if source_url is None:
        source_url = 'https://www.celestrak.com/NORAD/elements/'

    sats = {}
    _sats = load.tle(source_url + 'active.txt', reload=False)
    for xx in list(_sats.keys()): 
        if not isinstance(xx, int):
            sats[xx] = _sats[xx]
    return sats

def __get_sat_coord__(key, obs_time, obs_location, sats):

    coords = []
    #for key in keys:
    satellite = sats[key]
    topcentric = (satellite - obs_location).at(obs_time)
    alt, az, distance = topcentric.altaz()
    coords.append([az.radians, alt.radians])
    return [coords, key]

def get_sat_coods(sats, obs_time_list, obs_location):
    #for obs_time in obs_time_list:
    #    obs_time = obs_time.datetime.timetuple()
    #    obs_time = load.timescale().utc(*obs_time[:6])
    obs_time = obs_time_list[0]
    obs_time_range = [x.unix - obs_time.unix for x in obs_time_list]

    obs_time = obs_time.datetime.timetuple()
    obs_time = load.timescale(builtin=True).utc(*obs_time[:5] + (obs_time_range,))

    coords = []
    name = []

    pool =  Pool(16)
    r = pool.map(partial(__get_sat_coord__, sats=sats, 
            obs_time=obs_time, obs_location=obs_location), list(sats.keys()))
    for _r in r:
        coords.append(_r[0][0])
        name.append(_r[1])
    pool.close()
    pool.join()

    name = np.array(name)
    coords = np.array(coords)
    coords = np.rollaxis(coords, -1, 0)
    coords = np.ma.masked_invalid(coords)
    return coords, name

def remove_sats_below_horizen(sats, st_time, ed_time, location):

    #t1 = time.time()

    time_range = ed_time.unix - st_time.unix
    #time_range = st_time + np.linspace(0, time_range, 10) * u.second
    time_range = np.linspace(0, time_range, 10)
    obs_time = st_time.datetime.timetuple()
    obs_time = load.timescale(builtin=True).utc(*obs_time[:5] + (time_range, ))

    for key in list(sats.keys()):
        satellite = sats[key]
        topcentric = (satellite - location).at(obs_time)
        alt, az, distance = topcentric.altaz()
        if np.all(alt.radians < 0):
            del sats[key]

    #print "remove sat use %f s"%(time.time() - t1)
    return sats

class Satellite_Catalogue(object):


    def __init__(self, source_url=None, sats_type=None, reload=False):

        self.sats_type = sats_type
        if self.sats_type is None:
            self.sats_type =  _comms_tles + _gps_tles
        sats = []
        for sats_type in self.sats_type:
            sats.append(get_sat_tles(sats_type, reload=reload, source_url=source_url))
        self.sats = sats

        self._obs_time = None
        self._obs_time_list = np.zeros(1) * u.second
        self._obs_location = None

    @property
    def obs_time(self):
        return self._obs_time

    @obs_time.setter
    def obs_time(self, _obs_time):
        '''
        _obs_time UTC str
        '''
        if isinstance(_obs_time, str):
            self._obs_time = [Time(_obs_time, scale='utc'), ]
        elif isinstance(_obs_time, list):
            self._obs_time = [Time(_o, scale='utc') for _o in _obs_time]

    @property
    def obs_time_list(self):
        return self._obs_time_list

    @obs_time_list.setter
    def obs_time_list(self, _time_list):
        '''
        _time_list in astrou.units, ndarray or list of ndarray.
        if it is ndarray, same obs_time_list applied to all obs
        '''
        if self.obs_time is None:
            raise ValueError('Set obs_time firstly')

        if isinstance(_time_list, np.ndarray):
            self._obs_time_list = [_time_list, ] * len(self.obs_time)
        elif isinstance(_time_list, list):
            if len(_time_list) != len(self.obs_time):
                raise ValueError('time_list num does not match obs num')
            self._obs_time_list = _time_list

    @property
    def obs_location(self):
        return self._obs_location

    @obs_location.setter
    def obs_location(self, _location):
        '''
        _location = [Lat, Lon] in u.deg
        '''
        Lat, Lon = _location
        if np.sign(Lat.value) == -1:
            _Lat = "%f S"%abs(Lat.value)
        else:
            _Lat = "%f N"%abs(Lat.value)

        if np.sign(Lon.value) == -1:
            _Lon = "%f W"%abs(Lon.value)
        else:
            _Lon = "%f E"%abs(Lon.value)
        self._obs_location = Topos(_Lat, _Lon)

    def get_sate_coords(self):

        if self.obs_time is None:
            msg = 'Need to set obs_time'
            raise ValueError(msg)

        if self.obs_location is None:
            msg = 'Need to set obs_location'
            raise ValueError(msg)
    
        #obs_time_list = self.obs_time + self.obs_time_list
        coord_list_total = []
        name_list_total = []
        for oo, obs_start_time in enumerate(self.obs_time):
            obs_time_list = obs_start_time + self.obs_time_list[oo]
            coord_list = []
            name_list  = []

            print("Time range %s - %s"%(obs_time_list[0].utc, obs_time_list[-1].utc))
            for ii, sat_type in enumerate(self.sats_type):
                t0 = time.time()
                sats = copy.copy(self.sats[ii])
                #if sat_type == 'geo':
                #    sats = remove_sats_below_horizen(sats, obs_time_list[0], 
                #            obs_time_list[-1], self.obs_location)
                print("Satellite %12s has %4d satellites"%(sat_type, 
                                                         len(list(sats.keys()))), end=' ')
                
                coord_list_sats, name_list_sats\
                        = get_sat_coods(sats, obs_time_list, self.obs_location)
                coord_list.append(coord_list_sats)
                name_list.append(name_list_sats)
                print(" [use %6.2f s]"%(time.time() - t0))
            coord_list_total.append(coord_list)
            name_list_total.append(name_list)

        self.coord_list = coord_list_total
        self.name_list  = name_list_total

    def itersats_temperature_oneobs(self, obs_idx, pointings, beam_func=None):

        if isinstance(pointings, np.ndarray):
            'uisng same pointing for all obs'
            pointings = [pointings, ] * len(self.coord_list)
        elif isinstance(pointings, list):
            if len(pointings) != len(self.coord_list):
                raise ValueError('pointing num does not match obs num')

        oo = obs_idx

        pointing_Az  = pointings[oo][:, 0][:, None] * np.pi / 180.
        pointing_Alt = pointings[oo][:, 1][:, None] * np.pi / 180.

        for ii in range(len(self.sats_type)):

            coords = self.coord_list[oo][ii]
            if len(coords) == 0: continue
            coords = np.array(coords)
            coords_Az  = coords[:, :, 0]
            coords_Alt = coords[:, :, 1]

            sats_name  = self.name_list[oo][ii]

            _angle = np.sin(pointing_Alt) * np.sin(coords_Alt)\
                   + np.cos(pointing_Alt) * np.cos(coords_Alt)\
                   * np.cos(pointing_Az - coords_Az)
            _angle = np.arccos(_angle) * 180. / np.pi
            if beam_func is not None:
                _angle = beam_func(_angle)

            yield self.sats_type[ii], sats_name, _angle

    def itersats_temperature(self, pointings, beam_func=None):

        '''
        pointings : ndarray N x 2 (Az, Alt) in deg

        yield: separation angle if beam_func is None; or return the 
                beam convolved temperature, assuming the sats temperature
                is 1 K.
        '''
        if isinstance(pointings, np.ndarray):
            'uisng same pointing for all obs'
            pointings = [pointings, ] * len(self.coord_list)
        elif isinstance(pointings, list):
            if len(pointings) != len(self.coord_list):
                raise ValueError('pointing num does not match obs num')

        for oo, obs_start_time in enumerate(self.obs_time):

            obs_time_list = obs_start_time.unix\
                    + self.obs_time_list[oo].to(u.second).value
            obs_time_list = Time(obs_time_list, format='unix')
            obs_time_list.format = 'fits'

            for _info in self.itersats_temperature_oneobs(oo, pointings, beam_func):

                yield (obs_time_list,) + _info

    def get_angular_separation(self, pointings, beam_func=None):
        '''
        pointings : ndarray N x 2 (Az, Alt) in deg

        return: separation angle if beam_func is None; or return the 
                beam convolved temperature, assuming the sats temperature
                is 1 K.
        '''

        if isinstance(pointings, np.ndarray):
            'uisng same pointing for all obs'
            pointings = [pointings, ] * len(self.coord_list)
        elif isinstance(pointings, list):
            if len(pointings) != len(self.coord_list):
                raise ValueError('pointing num does not match obs num')

        self.angle_separation_list = []
        for oo, coord_list in enumerate(self.coord_list):

            pointing_Az  = pointings[oo][:, 0][:, None] * np.pi / 180.
            pointing_Alt = pointings[oo][:, 1][:, None] * np.pi / 180.

            angle_separation_list = []
            for ii in range(len(self.sats_type)):

                coords = coord_list[ii]
                if len(coords) == 0: continue
                coords = np.array(coords)
                coords_Az  = coords[:, :, 0]
                coords_Alt = coords[:, :, 1]

                _angle = np.sin(pointing_Alt) * np.sin(coords_Alt)\
                       + np.cos(pointing_Alt) * np.cos(coords_Alt)\
                       * np.cos(pointing_Az - coords_Az)
                _angle = np.arccos(_angle) * 180. / np.pi
                if beam_func is not None:
                    _angle = beam_func(_angle)
                angle_separation_list.append(_angle)

            self.angle_separation_list.append(angle_separation_list)

    def check_angular_separation(self, pointings, max_angle=10, beam_func=None,
            ymin=None, ymax=None, axes=None):

        #self.get_angular_separation(pointings, beam_func=beam_func)
        for oo, obs_start_time in enumerate(self.obs_time):
            if axes is None:
                fig = plt.figure(figsize=(8, 4))
                ax  = fig.add_axes([0.1, 0.1, 0.85, 0.85])
            else:
                ax = axes
                fig = ax.get_figure()

            obs_time_list = obs_start_time.unix\
                    + self.obs_time_list[oo].to(u.second).value
            x_axis = [ datetime.fromtimestamp(s) for s in obs_time_list]
            x_axis = mdates.date2num(x_axis)
            _angle_sum_total = np.zeros(len(obs_time_list))
            for _info in self.itersats_temperature_oneobs(oo, pointings, beam_func):
                sats_type, sats_name, _angle = _info
                ii = self.sats_type.index(sats_type)
                if len(_angle.shape) == 3:
                    # average over freq
                    _angle = np.mean(_angle, axis=0)

                ax.plot(x_axis, _angle, '-', color=_c_list[ii], lw=0.1)

                _angle_sum = np.sum(_angle, axis=1)
                _angle_sum_total += _angle_sum
                ax.plot(x_axis, _angle_sum, '-', color=_c_list[ii], ms=1, lw=1,
                        zorder=1000)

            ax.plot(x_axis, _angle_sum_total, '-', color='k', lw=1, zorder=2000)

            date_format = mdates.DateFormatter('%H:%M')
            ax.xaxis.set_major_formatter(date_format)
            ax.set_xlim(xmin=x_axis[0], xmax=x_axis[-1])
            ax.set_ylim(ymin=ymin, ymax=ymax)
            if beam_func is not None:
                ax.semilogy()
            ax.minorticks_on()
            ax.tick_params(length=4, width=1, direction='in')
            ax.tick_params(which='minor', length=2, width=1, direction='in')
        if axes is None:
            plt.show()

    def check_altaz(self):
        #coord_list = self.coord_list

        _n  = len(self.coord_list)
        fig = plt.figure(figsize=(4 * _n + 2 , 6))
        gs = gridspec.GridSpec(1, _n, right=0.95, wspace=0.1, figure=fig)
        #ax_r = fig2.add_subplot(gs[0,0])

        for oo, coord_list in enumerate(self.coord_list):
            #name_list = self.name_list[oo]
            #fig = plt.figure(figsize=(6, 6))
            #ax  = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='polar',
            #                    theta_offset=0.5 * np.pi, theta_direction=-1)
            ax = fig.add_subplot(gs[0,oo], projection='polar',
                    theta_offset=0.5 * np.pi, theta_direction=-1)
            
            legend_list = []
            for ii in range(len(self.sats_type)):

                
                coords = coord_list[ii]
                if len(coords) == 0: continue
                coords = np.array(coords)

                #names  = name_list[ii]

                for jj in range(coords.shape[1]):
                    good = coords[:, jj, 1] > 0
                    ax.plot(coords[good,jj,0], 90.- coords[good,jj,1] * 180./np.pi,
                            '.-', color=_c_list[ii], mfc=_c_list[ii], mec='none',
                            lw=0.1, ms=1)
                #ax.scatter(coords[:,0], 90.- coords[:,1] * 180./np.pi, s=10,
                #      marker='s', facecolor='none', edgecolor=_c_list[ii])

                legend_list.append(mpatches.Patch(color=_c_list[ii], 
                                                  label='%s'%self.sats_type[ii]))

            ax.set_rlim(0, 90)
            ax.set_title(self.obs_time[oo])
            if oo == 0:
                ax.legend(handles=legend_list, frameon=False, markerfirst=False, loc=1,
                    bbox_to_anchor=(1, 0.1, 0.15, 1))
        plt.show()


    def check_pointing(self, rising_altaz, setting_altaz, az_range=18., pad=2.):
        '''
        rising_altaz = [Az, Alt]
        setting_altaz = [Az, Alt]
        '''
        _az0, _alt0 = rising_altaz
        if _az0 < 0: _az0 += 360.

        _az1, _alt1 = setting_altaz
        if _az1 < 0: _az1 += 360.
            

        fig2 = plt.figure(figsize=(12, 12))
        gs = gridspec.GridSpec(1, 2, right=0.95, wspace=0.05, figure=fig2)
        ax_r = fig2.add_subplot(gs[0,0])
        ax_s = fig2.add_subplot(gs[0,1])

        coord_list_1 = self.coord_list[0]
        coord_list_2 = self.coord_list[1]
        name_list_1  = self.name_list[0]
        name_list_2  = self.name_list[1]

        _sats_list1 = []
        _sats_list2 = []

        legend_list = []
        for ii in range(len(self.sats_type)):
            coords_1 = coord_list_1[ii]
            if len(coords_1) == 0: continue
            coords_1 = np.array(coords_1)

            coords_2 = coord_list_2[ii]
            if len(coords_2) == 0: continue
            coords_2 = np.array(coords_2)

            names_1 = name_list_1[ii]
            names_2 = name_list_2[ii]

            for jj in range(coords_1.shape[1]):
                good  = coords_1[:, jj, 1] > 0
                az  = coords_1[good, jj, 0] * 180./np.pi
                alt = coords_1[good, jj, 1] * 180./np.pi
                #az[az>180.] -= 360.
                ax_r.plot(az, alt, '.-', color=_c_list[ii], mfc=_c_list[ii],
                    mec='none', lw=0.1, ) #label=names_1[jj])

                good  = coords_1[:, jj, 0]* 180./np.pi > _az0 - pad
                good *= coords_1[:, jj, 0]* 180./np.pi < _az0 + az_range + pad
                good *= coords_1[:, jj, 1]* 180./np.pi > _alt0 - pad - 0.5
                good *= coords_1[:, jj, 1]* 180./np.pi < _alt0 + pad + 0.5
                if np.any(good):
                    az  = np.mean(coords_1[good, jj, 0] * 180./np.pi)
                    alt = np.mean(coords_1[good, jj, 1] * 180./np.pi)# + shift1
                    _sats_list1.append([az, alt, r'%s'%names_1[jj], _c_list[ii]])

            for jj in range(coords_2.shape[1]):
                good = coords_2[:, jj,1] > 0
                az  = coords_2[good, jj, 0] * 180./np.pi
                alt = coords_2[good, jj, 1] * 180./np.pi
                ax_s.plot(az, alt, '.-', color=_c_list[ii], mfc=_c_list[ii],
                        mec='none', lw=0.1)

                good  = coords_2[:, jj, 0]* 180./np.pi > _az1 - pad
                good *= coords_2[:, jj, 0]* 180./np.pi < _az1 + az_range + pad
                good *= coords_2[:, jj, 1]* 180./np.pi > _alt1 - pad - 0.5
                good *= coords_2[:, jj, 1]* 180./np.pi < _alt1 + pad + 0.5
                if np.any(good):
                    az  = np.mean(coords_2[good, jj, 0] * 180./np.pi)
                    alt = np.mean(coords_2[good, jj, 1] * 180./np.pi) 
                    _sats_list2.append([az, alt, r'%s'%names_2[jj], _c_list[ii]])

            legend_list.append(mpatches.Patch(color=_c_list[ii], 
                label='%s'%self.sats_type[ii]))
        shift = 0.015
        _sats_list1.sort( key=lambda x: x[1])
        for ss in range(len(_sats_list1)):
            _sats = _sats_list1[ss]
            ax_r.annotate(_sats[2].replace('&', r'\&'), (_sats[0], _sats[1]), 
                    textcoords = 'axes fraction',
                    xytext = (0.8, 0.02 + ss * shift),
                    arrowprops=dict(edgecolor=_sats[3],
                        arrowstyle='->', 
                        connectionstyle="arc,angleA=180,armA=100,rad=20"),
                    horizontalalignment='left', 
                    verticalalignment='bottom',
                    size = 8,
                    color=_sats[3])

        _sats_list2.sort( key=lambda x: x[1])
        for ss in range(len(_sats_list2)):
            _sats = _sats_list2[ss]
            ax_s.annotate(_sats[2].replace('&', r'\&'), (_sats[0], _sats[1]), 
                    textcoords = 'axes fraction',
                    xytext = (0.8, 0.02 + ss * shift),
                    arrowprops=dict(edgecolor=_sats[3],
                        arrowstyle='->', 
                        connectionstyle="arc,angleA=180,armA=100,rad=20"),
                    horizontalalignment='left', 
                    verticalalignment='bottom',
                    size = 8,
                    color=_sats[3])

        ax_s.legend(handles=legend_list, frameon=False, markerfirst=False,
                loc=1, )

        ax_r.set_title(self.obs_time[0])
        ax_s.set_title(self.obs_time[1])

        ax_r.legend()
        
        xx0 = np.linspace(_az0, _az0 + az_range, 100)
        yy0 = np.ones_like(xx0) * _alt0
        ax_r.fill_between(xx0, yy0+0.5, yy0-0.5, facecolor='none', edgecolor='r')

        ax_r.set_xlim( _az0 - 10.,  _az0 + az_range + 10.)
        #ax_r.set_yticklabels([])
        ax_r.set_ylabel(r'${\rm Alt}\, [^{\circ}]$')
        ax_r.set_xlabel(r'${\rm AZ}\, [^{\circ}]$')
        ax_r.tick_params(which='both', direction='in')
        
        xx1 = np.linspace(_az1, _az1 + az_range, 100)
        yy1 = np.ones_like(xx1) * _alt1
        ax_s.fill_between(xx1, yy1+0.5, yy1-0.5, facecolor='none', edgecolor='b')

        ax_s.set_xlim( _az1 - 10.,  _az1 + az_range + 10.)
        ax_s.set_xlabel(r'${\rm AZ}\, [^{\circ}]$')
        #ax_s.set_ylabel(r'${\rm Alt}\, [^{\circ}]$')
        ax_s.set_yticklabels([])
        ax_s.tick_params(which='both', direction='in')
        
        _alt = (_alt0 + _alt1) * 0.5
        ax_r.set_ylim(_alt - 5., _alt + 5.)
        ax_s.set_ylim(_alt - 5., _alt + 5.)
        plt.show()


class MeerKATsite_Satellite_Catalogue(Satellite_Catalogue):

    #def __init__(self, reload=False, source_url=None):
    def __init__(self, *args, **kwargs):

        super(MeerKATsite_Satellite_Catalogue, self).__init__(*args, **kwargs)

        meerKAT_Lon =  (21. + 26./60. + 38.00/3600.) * u.deg #
        meerKAT_Lat = -(30. + 42./60. + 47.41/3600.) * u.deg #
        self.obs_location = [meerKAT_Lat, meerKAT_Lon]

    def get_angular_separation(self, pointings, beam_func=None):

        if beam_func is None:
            # using modified sinc function, quite close to Khan's model
            beam_func = lambda x: np.sinc(x) ** 2 * ((1/(np.abs(x) + 1)) ** 2.5)

        super(MeerKATsite_Satellite_Catalogue, self).get_angular_separation(
                pointings, beam_func=beam_func)

    def check_angular_separation(self, pointings, max_angle=100, beam_func=None,
            ymin=None, ymax=None, axes=None):


        if beam_func is None:
            # using modified sinc function, quite close to Khan's model
            beam_func = lambda x: np.sinc(x) ** 2 * ((1/(np.abs(x) + 1)) ** 2.5)


        super(MeerKATsite_Satellite_Catalogue, self).check_angular_separation(
                pointings, max_angle=max_angle,
                beam_func=beam_func, ymin=ymin, ymax=ymax, axes=axes)



    def check_pointing(self, altaz, figaxes=None, pad=2.):
        '''
        rising_altaz = [Az, Alt]
        setting_altaz = [Az, Alt]
        '''
        #_az0, _alt0 = altaz
        _az_min, _alt_min = np.min(altaz, axis=0)
        _az_max, _alt_max = np.max(altaz, axis=0)
        if _az_min < 0: _az_min += 360.
        if _az_max < 0: _az_max += 360.

        if figaxes is None:
            fig2 = plt.figure(figsize=(12, 12))
            gs = gridspec.GridSpec(1, 1, right=0.95, wspace=0.05, figure=fig2)
            ax_r = fig2.add_subplot(gs[0,0])
        else:
            fig2, ax_r = figaxes

        coord_list_1 = self.coord_list[0]
        name_list_1  = self.name_list[0]

        _sats_list1 = []

        legend_list = []
        for ii in range(len(self.sats_type)):
            coords_1 = coord_list_1[ii]
            if len(coords_1) == 0: continue
            coords_1 = np.ma.masked_invalid(coords_1)

            names_1 = name_list_1[ii]

            for jj in range(coords_1.shape[1]):
                good  = coords_1[:, jj, 1] > 0
                az  = coords_1[good, jj, 0] * 180./np.pi
                alt = coords_1[good, jj, 1] * 180./np.pi
                #az[az>180.] -= 360.
                az[az<0] += 360.
                ax_r.plot(az, alt, '.-', color=_c_list[ii], mfc=_c_list[ii],
                    mec='none', lw=0.1, ) #label=names_1[jj])

                good  = coords_1[:, jj, 0]* 180./np.pi > _az_min - pad
                good *= coords_1[:, jj, 0]* 180./np.pi < _az_max + pad
                good *= coords_1[:, jj, 1]* 180./np.pi > _alt_min - pad - 0.5
                good *= coords_1[:, jj, 1]* 180./np.pi < _alt_max + pad + 0.5
                if np.any(good):
                    az  = np.mean(coords_1[good, jj, 0] * 180./np.pi)
                    alt = np.mean(coords_1[good, jj, 1] * 180./np.pi)
                    _sats_list1.append([az, alt, r'%s'%names_1[jj], _c_list[ii]])

            legend_list.append(mpatches.Patch(color=_c_list[ii], 
                label='%s'%self.sats_type[ii]))
        #shift = 0.015
        shift = 0.02
        _sats_list1.sort( key=lambda x: x[1])
        for ss in range(len(_sats_list1)):
            _sats = _sats_list1[ss]
            ax_r.annotate(_sats[2].replace('&', r'\&'), (_sats[0], _sats[1]), 
                    textcoords = 'axes fraction',
                    xytext = (0.8, 0.02 + ss * shift),
                    arrowprops=dict(edgecolor=_sats[3],
                        arrowstyle='->', 
                        connectionstyle="arc,angleA=180,armA=100,rad=20"),
                    horizontalalignment='left', 
                    verticalalignment='bottom',
                    size = 8,
                    color=_sats[3])

        ax_r.legend(handles=legend_list, frameon=False, markerfirst=False,
                loc=1, )

        ax_r.set_title(self.obs_time[0])

        xx0 = np.linspace(_az_min, _az_max, 100)
        yy0 = np.ones_like(xx0) * (_alt_min + _alt_max) * 0.5
        ax_r.fill_between(xx0, yy0+0.5, yy0-0.5, facecolor='none', edgecolor='r')

        ax_r.set_xlim( _az_min - 10.,  _az_max + 10.)
        ax_r.set_ylabel(r'${\rm Alt}\, [^{\circ}]$')
        ax_r.set_xlabel(r'${\rm AZ}\, [^{\circ}]$')
        ax_r.tick_params(which='both', direction='in')

        _alt = (_alt_min + _alt_max) * 0.5
        ax_r.set_ylim(_alt - 5., _alt + 5.)
        if figaxes is None:
            plt.show()
            plt.close()
