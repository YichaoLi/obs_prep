#!/usr/bin/env python
from pylab import rcParams
#rcParams['figure.facecolor'] = 'white'
#rcParams['figure.facecolor'] = 'white'
#rcParams['figure.dpi'] = 80

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.dates as mdates
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz
from astropy.coordinates import get_sun, get_moon, FK5
from astropy.coordinates import EarthLocation
from astroplan import Observer

import obs_prep.check_satellite as cs

import pandas as pd
import numpy as np
import copy

from IPython.display import HTML

_c_list = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
          "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]

class Calibrator(object):

    def __init__(self, ra, dec, flux=None, name='Cal', cal_time=0):

        self.ra  = ra
        self.dec = dec
        self.flux = flux
        self.name = name
        self.cal_time = cal_time

class SkyField(object):

    def __init__(self, ra_min, ra_max, dec_min, dec_max, name='SkyField'):

        self.ra_min = ra_min * u.deg
        self.ra_max = ra_max * u.deg
        self.dec_min = dec_min * u.deg
        self.dec_max = dec_max * u.deg
        self.name = name
        super(SkyField, self).__init__()

    @property
    def ra(self):
        return (self.ra_min + self.ra_max) * 0.5 

    @property
    def dec(self):
        return (self.dec_min + self.dec_max) * 0.5 

    @property
    def edges(self):
        ra_min = self.ra_min.to(u.deg).value
        ra_max = self.ra_max.to(u.deg).value
        dec_min = self.dec_min.to(u.deg).value
        dec_max = self.dec_max.to(u.deg).value

        edg = [ np.concatenate([
                np.ones(100) * ra_min, 
                np.linspace(ra_min, ra_max, 100), 
                np.ones(100) * ra_max,
                np.linspace(ra_max, ra_min, 100) ]) * u.deg,
                np.concatenate([
                np.linspace(dec_min, dec_max, 100),
                np.ones(100) * dec_max,
                np.linspace(dec_max, dec_min, 100),
                np.ones(100) * dec_min,]) * u.deg ]
        return edg

    def check_field(self, axes=None, project=False, c='k'):

        if axes is None:
            fig = plt.figure(figsize=(10, 6))
            ax  = fig.add_axes([0.1, 0.1, 0.85, 0.85])
        else:
            fig, ax = axes

        #print self.ra_min
        #print self.edges[0]
        #print self.edges[1]

        edg_ra  = self.edges[0].value
        edg_dec = self.edges[1].value
        if project:
            edg_ra  -= self.ra.to(u.deg).value
            edg_ra  *= np.cos(edg_dec * np.pi / 180.)
            edg_dec -= self.dec.to(u.deg).value

        ax.plot(edg_ra, edg_dec, c + '-', linewidth=1.5, zorder=100)

        ax.set_aspect('equal')
        if axes is None:
            ax.set_xlabel('RA')
            ax.set_ylabel('Dec')

        return fig, ax

class ObservationSchedule(object):

    def __init__(self):

        super(ObservationSchedule, self).__init__()

        self._obs_date_list = None
        self._location = None
        self._msg = ''
        self._obs_info = []

        self._observed_radec = []

    @property
    def obs_details(self):
        return self._msg

    @property
    def obs_info(self):
        return self._obs_info

    @property
    def obs_table_string(self):
        return[pd.DataFrame(_obs_info).to_string(index=False)
                for _obs_info in self._obs_info]

    @property
    def obs_table(self):
        return[HTML(pd.DataFrame(_obs_info).to_html(index=False))
                for _obs_info in self._obs_info]

    @property
    def obs_location(self):
        return self._location

    @obs_location.setter
    def obs_location(self, location):
        '''
        location: [Lon, Lat] in deg
        '''
        _Lon, _Lat = location
        self._location = EarthLocation.from_geodetic(_Lon, _Lat) 

    @property
    def obs_date(self):
        '''
        The date(s) for observation
        '''
        return self._obs_date_list

    @obs_date.setter
    def obs_date(self, date):
        '''
        date: the utc date time you are planing to obs
        '''

        if not isinstance(date, list):
            date = [date, ]

        _loc = self.obs_location
        if _loc is None:
            msg = 'location not specified'
            Warning(msg)

        self._obs_date_list = [Time(x, location=_loc) for x in date]

    @property
    def observed_radec(self):

        return self._observed_radec

    @observed_radec.setter
    def observed_radec(self, radec):

        if radec is None:
            self._observed_radec = []
        else:
            self._observed_radec.append(radec)



class ScaningStrategy(object):

    def __init__(self):

        super(ScaningStrategy, self).__init__()

        self._target_field_list = None

        self.dumping_rate = None

        self._obs_date_list = None
        self._block_time = None

        self.alt_list = None
        self.az_list  = None
        self.t_list   = None
        self.ra_list  = None
        self.dec_list = None
        self._cal_list = None

    def align_with_obsdate(self, value):

        if not isinstance(value, list) or not isinstance(value[0], list):
            value = [value, ] * len(self._obs_date_list)
        return value

    @property
    def target_field(self):
        return self._target_field_list

    @target_field.setter
    def target_field(self, field):
        if not isinstance(field, list):
            field = [field, ]
        self._target_field_list = field

    @property
    def calibrator_list(self):
        return self._cal_list

    @calibrator_list.setter
    def calibrator_list(self, cal):
        if not isinstance(cal, list):
            cal = [cal, cal, cal, cal]
        self._cal_list = self.align_with_obsdate(cal)

    def generate_altaz(self, location=None):
        pass

    def get_schedule_field(self,):

        pass

    def check_schedule(self):

        _length   = self.check_length.to(u.hour).value
        _step     = self.check_step.to(u.hour).value
        _location = self.obs_location
        #_alt_list    = self.obs_alt
        #_target_list = self.target_field
        #obs_az_range = self.az_raster_range 
        block_time   = self.block_time

        ndays = len(self.obs_date)
        fig = plt.figure(figsize=(8, 2 * ndays + 1))
        gs = gridspec.GridSpec(ndays, 1,
                left=0.12, bottom=0.08, right=0.9, top=0.95, wspace=0.1, hspace=0.0)

        for oo, obs_date in enumerate(self.obs_date):

            try:
                cal_list = self.calibrator_list
                if cal_list is not None:
                    cal_list = list(set(cal_list[oo]))
            except AttributeError:
                cal_list = None

            try: 
                alt_list = self.obs_alt
                if alt_list is not None:
                    alt_list = alt_list[oo]
            except AttributeError:
                alt_list = None

            try:
                alt_list = self.alt_limit
                if alt_list is not None:
                    alt_list = alt_list[oo]
            except AttributeError:
                alt_list = None


            ax = fig.add_subplot(gs[oo, 0])

            _obs_time = obs_date + np.arange(0, _length, _step) * u.hour
            field = self.target_field
            plot_sch_oneday(_obs_time, _location, field, ax, alt_list, cal_list)
            
            if oo == 0:
                ax.legend(frameon=False, ncol=3)
            if oo == ndays - 1:
                ax.set_xlabel('UTC')
            else:
                ax.set_xticklabels([])
        plt.show()

class DriftScan(ScaningStrategy):

    def __init__(self):

        super(ScaningStrategy, self).__init__()

        self._scaning_speed = None

        self._az_raster_range = None
        self._obs_alt = None

        self.start_time_list = None
        self.start_az_list   = None
        self.start_alt_list  = None

    @property
    def scaning_speed(self):
        return self._scaning_speed
    
    @scaning_speed.setter
    def scaning_speed(self, value):
        if not isinstance(value, list):
            value = [value, value]
        self._scaning_speed = self.align_with_obsdate(value)

    @property
    def obs_alt(self):
        return self._obs_alt

    @obs_alt.setter
    def obs_alt(self, alt):
        if not isinstance(alt, list):
            alt = [alt, alt]
        self._obs_alt = self.align_with_obsdate(alt)

    @property
    def block_time(self):
        return self._block_time

    @block_time.setter
    def block_time(self, blocktime):
        if not isinstance(blocktime, list):
            blocktime = [blocktime, blocktime]
        self._block_time = self.align_with_obsdate(blocktime)

    def generate_altaz(self, location=None, shift=0):

        pass

    def check_scans(self, axes=None, project=False, day_select=None, 
            ra_range=[None, None], dec_range=[None, None], zero_center=False,
            combine_plot=True):

        _location = self.obs_location
        obs_date = np.arange(len(self.obs_date))
        if day_select is not None:
            obs_date = obs_date[slice(*day_select)]

        _c = ['r', 'b']

        for ff, target_field in enumerate(self.target_field):
            fig, ax = target_field.check_field(axes=axes, project=project)
            if combine_plot: axes = (fig, ax)
            dd = 0
            for ii in obs_date:
                for jj in range(2):
                    _alt = self.alt_list[ff][ii][jj]
                    _az  = self.az_list[ff][ii][jj]
                    _t   = self.t_list[ff][ii][jj]

                    pp = SkyCoord(alt=_alt, az=_az, frame='altaz', 
                                  location=_location, obstime=_t)
                    pp = pp.transform_to('icrs')

                    ra  = pp.ra.deg
                    dec = pp.dec.deg
                    if zero_center:
                        ra[ra>180] -= 360.

                    if project:
                        ra -= target_field.ra.value
                        ra *= np.cos(target_field.dec.to(u.rad).value)
                        dec -= target_field.dec.value

                    ax.plot(ra, dec, _c[jj]+'.-', ms=1., mfc='none', linewidth=0.5)
                    ax.plot(ra[0], dec[0], 'go')
                    ax.text(ra[0], dec[0] - 0.2, '%d'%dd)
                    dd += 1
            for _radec in self.observed_radec:
                ra  = _radec[:, 0].copy()
                dec = _radec[:, 1].copy()
                if project:
                    ra -= target_field.ra.value
                    ra *= np.cos(target_field.dec.to(u.rad).value)
                    dec -= target_field.dec.value
                ax.plot(ra, dec, 'g.-', ms=1., mfc='none', 
                        linewidth=0.5, zorder=10000)
            ax.set_xlim(tuple(ra_range))
            ax.set_ylim(tuple(dec_range))
            if axes is None:
                plt.show()

    def check_hitmap(self, axes=None, day_select=None, pad=(3, 1), 
            zero_center=False, combine_plot=True, edges=None):

        project=False

        _location = self.obs_location
        obs_date = np.arange(len(self.obs_date))
        if day_select is not None:
            obs_date = obs_date[slice(*day_select)]

        _c = ['r', 'b']

        hitmap = None
        for ff, target_field in enumerate(self.target_field):
            fig, ax = target_field.check_field(axes=axes, project=project)
            if combine_plot: axes = (fig, ax)
            pixsize = 0.3
            if edges is None:
                ra_min  = target_field.ra_min.value
                ra_max  = target_field.ra_max.value
                dec_min = target_field.dec_min.value
                dec_max = target_field.dec_max.value
            else:
                ra_min, ra_max, dec_min, dec_max = edges

            ra_range  = [ra_min, ra_max]
            dec_range = [dec_min, dec_max]

            ra_bins_e  = np.arange(ra_min - pad[0] - 0.5*pixsize, 
                                   ra_max + pad[0] + pixsize, 
                                   pixsize)
            dec_bins_e = np.arange(dec_min- pad[1] - 0.5*pixsize, 
                                   dec_max+ pad[1] + pixsize,
                                   pixsize)
            ra_bins  = ra_bins_e[:-1] + 0.5*pixsize
            dec_bins = dec_bins_e[:-1] + 0.5*pixsize
            #if (not combine_plot) or (hitmap is None):
            if hitmap is None:
                hitmap = np.zeros(ra_bins.shape + dec_bins.shape)
            dd = 0
            for ii in obs_date:
                for jj in range(2):
                    _alt = self.alt_list[ff][ii][jj]
                    _az  = self.az_list[ff][ii][jj]
                    _t   = self.t_list[ff][ii][jj]

                    pp = SkyCoord(alt=_alt, az=_az, frame='altaz', 
                                  location=_location, obstime=_t)
                    pp = pp.transform_to('icrs')

                    ra  = pp.ra.deg
                    dec = pp.dec.deg
                    if zero_center:
                        ra[ra>180] -= 360.

                    hitmap += np.histogram2d(ra, dec, bins=[ra_bins_e, dec_bins_e])[0]

                    dd += 1

            if not combine_plot:
                hitmap = np.ma.masked_equal(hitmap, 0)
                vmin = 0
                vmax = hitmap.max() * 1.5
                ax.pcolormesh(ra_bins_e, dec_bins_e, hitmap.T, 
                        vmin=vmin, vmax=vmax, cmap='Blues')

                ax.set_xlim(tuple(ra_range))
                ax.set_ylim(tuple(dec_range))
            if axes is None:
                plt.show()

        if combine_plot:
            hitmap = np.ma.masked_equal(hitmap, 0)
            vmin = 0
            vmax = hitmap.max() * 1.1
            ax.pcolormesh(ra_bins_e, dec_bins_e, hitmap.T, 
                    vmin=vmin, vmax=vmax, cmap='Blues')

            ax.set_xlim(tuple(ra_range))
            ax.set_ylim(tuple(dec_range))

class HorizonRasterDrift(DriftScan):

    @property
    def az_raster_range(self):
        return self._az_raster_range

    @az_raster_range.setter
    def az_raster_range(self, az_range):
        if not isinstance(az_range, list):
            az_range = [az_range, az_range]
        self._az_raster_range = self.align_with_obsdate(az_range)

    def generate_altaz(self, location=None, shift=0.):

        obs_speed    = self.scaning_speed 
        obs_int      = self.dumping_rate 
        block_time   = self.block_time
        obs_az_range = self.az_raster_range #self.params['obs_az_range']

        self.alt_list = []
        self.az_list  = []
        self.t_list   = []
        start_time_list = self.start_time_list
        start_az_list  = self.start_az_list
        start_alt_list = self.start_alt_list

        nfields = len(self.target_field)
        ndays = len(self.obs_date)
        #t_delay_pd = (obs_az_range[0][0] * 2. / obs_speed).decompose() / float(ndays)
        #t_delay_pd = (np.arange(ndays) - 0.5 * (ndays - 1)) * t_delay_pd


        for ff in range(nfields):
            alt_list = []
            az_list  = []
            t_list   = []
            for ii in range(ndays):
                day_alt_list = []
                day_az_list  = []
                day_t_list   = []
                for jj in range(2):
                    _obs_speed = obs_speed[ii][jj]
                    t_delay_pd  = (obs_az_range[ii][jj] * 2. / _obs_speed).decompose()
                    t_delay_pd /= float(ndays)
                    t_delay_pd  = (ii - 0.5 * (ndays - 1)) * t_delay_pd
                    starttime = start_time_list[ff][ii][jj] + t_delay_pd
                    alt_start = start_alt_list[ff][ii][jj]
                    az_start  = start_az_list[ff][ii][jj]

                    obs_len   = int((block_time[ii][jj] / obs_int).decompose().value)

                    _alt_list = (np.ones(obs_len) * alt_start).value

                    day_alt_list.append(np.array(_alt_list) * u.deg)

                    _az_space = (_obs_speed * obs_int).to(u.deg)
                    _one_way_npoints = \
                            (obs_az_range[ii][jj] / _obs_speed / obs_int).decompose()
                    _az_list = np.arange(_one_way_npoints)
                    _az_list = np.append(_az_list, _az_list[::-1] + 1)
                    _az_list = _az_list * _az_space
                    _az_list += az_start
                    if jj == 0: _az_list -= shift * u.deg
                    elif jj == 1: _az_list += shift * u.deg
                    _az_list = _az_list.value
                    _az_list = [_az_list[i%int(2.*_one_way_npoints)] 
                            for i in range(obs_len)]

                    day_az_list.append(np.array(_az_list) * u.deg)

                    _time = np.arange(obs_len) * obs_int + starttime
                    _slew = (np.round(np.arange(obs_len)/_one_way_npoints))*5.*u.second
                    _time += _slew
                    day_t_list.append(Time(_time, location=location))

                alt_list.append(day_alt_list)
                az_list.append(day_az_list)
                t_list.append(day_t_list)

            self.alt_list.append(alt_list)
            self.az_list.append(az_list)
            self.t_list.append(t_list)

        #self.alt_list = np.concatenate(alt_list) * u.deg
        #self.az_list  = np.concatenate(az_list) * u.deg
        #self.t_list   = Time(np.concatenate(t_list))
        #print self.alt_list
        #self.alt_list = np.array(self.alt_list) * u.deg
        #self.az_list  = np.array(self.az_list) * u.deg


class MeerKATHRD(ObservationSchedule, HorizonRasterDrift):

    def __init__(self):

        super(MeerKATHRD, self).__init__()

        self.check_step = 0.1 * u.hour
        self.check_length = 24 * u.hour

        meerKAT_Lon =  (21. + 26./60. + 38.00/3600.) * u.deg #
        meerKAT_Lat = -(30. + 42./60. + 47.41/3600.) * u.deg #
        self.obs_location = [meerKAT_Lon, meerKAT_Lat]

    def generate_altaz(self, shift=0):

        super(MeerKATHRD, self).generate_altaz(self.obs_location, shift=shift)

    def get_schedule(self):

        self.start_time_list = []
        self.start_az_list   = []
        self.start_alt_list  = []
        for target_field in self.target_field:

            field_time_list, field_az_list, field_alt_list = \
                    self.get_schedule_field(target_field, center_shift = (0.53, 0.47))
            
            self.start_time_list.append(field_time_list)
            self.start_az_list.append(field_az_list)
            self.start_alt_list.append(field_alt_list)



    def get_schedule_field(self, target_field, center_shift = (0.53, 0.47)):
        

        _length   = self.check_length.to(u.hour).value
        _step     = self.check_step.to(u.hour).value
        _location = self.obs_location
        _alt_list    = self.obs_alt
        _target_list = target_field
        obs_az_range = self.az_raster_range 
        block_time   = self.block_time

        field_start_time_list = []
        field_start_az_list   = []
        field_start_alt_list  = []
        for oo,  obs_date in enumerate(self.obs_date):

            _alt = _alt_list[oo]
            _target = SkyCoord(ra = _target_list.ra, dec = _target_list.dec)

            _obs_time = obs_date + np.arange(0, _length, _step) * u.hour

            rt, st = _get_risingtime_at_alt(_obs_time, _location, _alt, _target)

            _obs_time = rt + (np.linspace(-1, 1, 100) * _step) * u.hour
            rt = _get_risingtime_at_alt(_obs_time, _location, _alt, _target)[0]
            _raltaz = _target.transform_to(AltAz(obstime=rt, location=_location))
            ralt, raz = _raltaz.alt.deg, _raltaz.az.deg
            raz -= center_shift[0] * obs_az_range[oo][0].to(u.deg).value
            # make sure the az range between -185 to 275
            if raz < -185.: raz += 360.
            if raz + obs_az_range[oo][0].to(u.deg).value > 275. :  raz -= 360.
            rt = rt - block_time[oo][0].to(u.hour) * 0.5

            _obs_time = st + (np.linspace(-1, 1, 100) * _step) * u.hour
            st = _get_risingtime_at_alt(_obs_time, _location, _alt, _target)[1]
            _saltaz = _target.transform_to(AltAz(obstime=st, location=_location))
            salt, saz = _saltaz.alt.deg, _saltaz.az.deg
            saz -= center_shift[1] * obs_az_range[oo][1].to(u.deg).value
            # make sure the az range between -185 to 275
            if saz < -185.: saz += 360.
            if saz + obs_az_range[oo][1].to(u.deg).value > 275. :  saz -= 360.
            st = st - block_time[oo][1].to(u.hour) * 0.5

            field_start_time_list.append([rt, st])
            field_start_az_list.append(  [raz  * u.deg, saz  * u.deg])
            field_start_alt_list.append( [ralt * u.deg, salt * u.deg])

        return field_start_time_list, field_start_az_list, field_start_alt_list


    def check_schedule(self):

        _length   = self.check_length.to(u.hour).value
        _step     = self.check_step.to(u.hour).value
        _location = self.obs_location
        _alt_list    = self.obs_alt
        #_target_list = self.target_field
        obs_az_range = self.az_raster_range 
        block_time   = self.block_time

        ndays = len(self.obs_date)
        fig = plt.figure(figsize=(8, 2 * ndays + 1))
        gs = gridspec.GridSpec(ndays, 1,
                left=0.12, bottom=0.08, right=0.9, top=0.95, wspace=0.1, hspace=0.0)

        for oo, obs_date in enumerate(self.obs_date):

            cal_list = list(set(self.calibrator_list[oo]))

            ax = fig.add_subplot(gs[oo, 0])

            _obs_time = obs_date + np.arange(0, _length, _step) * u.hour
            field = self.target_field
            plot_sch_oneday(_obs_time, _location, field, ax, _alt_list[oo], cal_list)

            if oo == 0:
                ax.legend(frameon=False, ncol=3)
            if oo == ndays - 1:
                ax.set_xlabel('UTC')
        plt.show()

    def make_obs(self, output=None):

        _location = self.obs_location
        obs_speed = self.scaning_speed 
        obs_int   = self.dumping_rate 
        obs_az_range = self.az_raster_range 
        block_time   = self.block_time

        self._msg = ''

        self._msg += '='*60
        self._msg += '\n\n'
    
        self._msg +=  "# Log.Lat.   %12.7f %12.7f\n"%(_location.lon.deg,_location.lat.deg)
        self._msg +=  "# Int.Time   %10.5f s\n"%(obs_int.to(u.second).value)
        self._msg +=  '-'*60
        self._msg +=  '\n\n'

        msg = ''
        msg +=  "Log.Lat.   %12.7f %12.7f\n"%(_location.lon.deg,_location.lat.deg)
        msg +=  "AZRange    %f deg\n"%(obs_az_range[0][0].value)
        msg +=  "ScanSpeed  %f arcmin / s\n"%(obs_speed[0][0].value)
        msg +=  "Int.Time   %10.5f s\n"%(obs_int.to(u.second).value)
        msg +=  "SlewTime   0.000 s\n"
        msg +=  '-'*60

        #for key in wgz_field.keys():
            #fig = plt.figure(figsize=(10, 4))
            #ax  = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        self._obs_info = []
        for ff, field in enumerate( self.target_field ):
            _obs_info, _msg = make_obs_oneday(obs_az_range, obs_speed, self.obs_date,
                    field, self.alt_list[ff], self.az_list[ff], self.t_list[ff],
                    self.calibrator_list, block_time, _location)
            self._obs_info.append(_obs_info)
            self._msg += _msg

            _info = []
            for _l in _obs_info:
                az, alt = _l[9].split()
                _info.append((_l[0][0], _l[7], _l[8], float(az), float(alt), _l[3] ))
            _info = np.array(_info, dtype=[
                ('Obs.', 'S2'), 
                ('Start Time [UTC]', 'S25'),
                ('Start Time [LST]', 'f4'),
                ('Pointing (AZ)', 'f4'),
                ('Pointing (Alt)', 'f4'),
                ('Total Time [s]', 'f4'),])
            if output is not None:
                np.savetxt(output + '_%s.dat'%field.name, _info,
                           fmt='%s,%s,%12.7f,%12.7f,%12.7f,%f', delimiter=',',
                           header=msg)

    def check_satellite(self, reload=False, space=10):

        msc = cs.MeerKATsite_Satellite_Catalogue(reload=reload)
        for ff, field in enumerate( self.target_field ):
            for oo, obs_date in enumerate(self.obs_date):
                rt = self.t_list[ff][oo][0][::space]
                st = self.t_list[ff][oo][1][::space]
                msc.obs_time = [rt[0], st[0]]
                msc.obs_time_list = [rt - rt[0], st - st[0]]
                msc.get_sate_coords()

                r_az  = self.az_list[ff][oo][0][::space].to(u.deg).value
                r_alt = self.alt_list[ff][oo][0][::space].to(u.deg).value

                s_az  = self.az_list[ff][oo][1][::space].to(u.deg).value
                s_alt = self.alt_list[ff][oo][1][::space].to(u.deg).value

                msc.check_altaz()
                #msc.check_pointing([r_az[0], r_alt[0]], [s_az[0], s_alt[0]],
                #        az_range = self.az_raster_range[oo][0].to(u.deg).value)

                fig2 = plt.figure(figsize=(12, 12))
                gs = gridspec.GridSpec(1, 2, right=0.95, wspace=0.05, figure=fig2)
                ax_r = fig2.add_subplot(gs[0,0])
                ax_s = fig2.add_subplot(gs[0,1])

                msc.check_pointing(
                        np.concatenate([r_az[:, None], r_alt[:, None]], axis=1), 
                        figaxes=(fig2, ax_r))
                        #az_range = self.az_raster_range[oo][0].to(u.deg).value)
                msc.check_pointing(
                        np.concatenate([s_az[:, None], s_alt[:, None]], axis=1),
                        figaxes=(fig2, ax_s))
                        #az_range = self.az_raster_range[oo][0].to(u.deg).value)
                ax_s.set_yticklabels([])
                pointing = [np.concatenate([r_az[:, None], r_alt[:, None]], axis=1),
                            np.concatenate([s_az[:, None], s_alt[:, None]], axis=1),]
                msc.check_angular_separation(pointing, ymin=1.e-10, ymax=1)

class MultiDriftScan(DriftScan):

    def generate_altaz(self, location=None):

        obs_int      = self.dumping_rate
        block_time   = self.block_time

        self.alt_list = []
        self.az_list  = []
        self.t_list   = []
        start_time_list = self.start_time_list
        start_az_list  = self.start_az_list
        start_alt_list = self.start_alt_list

        nfields = len(self.target_field)
        ndays = len(self.obs_date)

        for ff in range(nfields):
            alt_list = []
            az_list  = []
            t_list   = []
            for ii in range(ndays):
                day_alt_list = []
                day_az_list  = []
                day_t_list   = []
                for jj in range(2):
                    #_obs_speed = obs_speed[ii][jj]
                    #t_delay_pd  = (obs_az_range[ii][jj] * 2. / _obs_speed).decompose()
                    #t_delay_pd /= float(ndays)
                    #t_delay_pd  = (ii - 0.5 * (ndays - 1)) * t_delay_pd
                    starttime = start_time_list[ff][ii][jj]
                    alt_start = start_alt_list[ff][ii][jj]
                    az_start  = start_az_list[ff][ii][jj]

                    obs_len   = int((block_time[ii][jj] / obs_int).decompose().value)

                    _alt_list = (np.ones(obs_len) * alt_start).value
                    day_alt_list.append(np.array(_alt_list) * u.deg)
                    _az_list = (np.ones(obs_len) * az_start).value
                    day_az_list.append(np.array(_az_list) * u.deg)

                    _time = np.arange(obs_len) * obs_int + starttime
                    #_slew = (np.round(np.arange(obs_len)/_one_way_npoints))*5.*u.second
                    #_time += _slew
                    day_t_list.append(Time(_time, location=location))

                alt_list.append(day_alt_list)
                az_list.append(day_az_list)
                t_list.append(day_t_list)

            self.alt_list.append(alt_list)
            self.az_list.append(az_list)
            self.t_list.append(t_list)

    def make_obs(self, output=None):

        _location = self.obs_location
        obs_int   = self.dumping_rate
        block_time   = self.block_time

        self._msg = ''

        #self._msg += '='*60
        #self._msg += '\n\n'

        self._msg +=  "Log.Lat.   %12.7f %12.7f\n"%(_location.lon.deg,_location.lat.deg)
        self._msg +=  "AZRange    0.000 deg\n"
        self._msg +=  "ScanSpeed  0.000 arcmin / s\n"
        self._msg +=  "Int.Time   %10.5f s\n"%(obs_int.to(u.second).value)
        self._msg +=  "SlewTime   0.000 s\n"
        self._msg +=  '-'*60
        #self._msg +=  '\n\n'

        self._obs_info = []


        #status = ['Rising Block', 'Setting Block']
        status = ['D', 'D']


        for ff, field in enumerate( self.target_field ):
            _obs_info = []
            for oo, obs_date in enumerate(self.obs_date):
                for ii in range(2):

                    _obs_tot  = block_time[oo][ii].to(u.second).value

                    start_time = self.t_list[ff][oo][ii][0]
                    end_time = self.t_list[ff][oo][ii][-1]
                    _alt_0 = self.alt_list[ff][oo][ii][0].to(u.deg).value
                    _az_0 = self.az_list[ff][oo][ii][0].to(u.deg).value

                    _obs_info.append((status[ii],
                        #'%s + %s'%(self.calibrator_list[oo][ii * 2].name,
                        #           self.calibrator_list[oo][ii * 2 + 1].name),
                        '%s'%start_time, '%s'%end_time,
                        _az_0, _alt_0, _obs_tot, ))

            _obs_info = np.array(_obs_info, dtype = [
                ('Obs.', 'S2'), #('Calibrator', 'S15'),
                ('Start Time [UTC]', 'S25'),
                ('End Time [UTC]', 'S25'),
                ('Pointing (AZ)', 'f4'),
                ('Pointing (Alt)', 'f4'),
                ('Total Time [s]', 'f4'),])
            self._obs_info.append(_obs_info)

            if output is not None:
                np.savetxt(output + '_%s.dat'%field.name, _obs_info,
                           fmt='%s,%s,%s,%6.2f,%6.2f,%f', delimiter=',',
                           header=self._msg)
class MeerKATMDS(ObservationSchedule, MultiDriftScan):

    def __init__(self):

        super(MeerKATMDS, self).__init__()

        self.check_step = 0.1 * u.hour
        self.check_length = 24 * u.hour

        meerKAT_Lon =  (21. + 26./60. + 38.00/3600.) * u.deg #
        meerKAT_Lat = -(30. + 42./60. + 47.41/3600.) * u.deg #
        self.obs_location = [meerKAT_Lon, meerKAT_Lat]

    def get_schedule(self, repeat=1):

        self.start_time_list = []

        self.start_az_list   = []
        self.start_alt_list  = []
        for target_field in self.target_field:

            field_time_list, field_az_list, field_alt_list = \
                    self.get_schedule_field(target_field, repeat=repeat)

            self.start_time_list.append(field_time_list)
            self.start_az_list.append(field_az_list)
            self.start_alt_list.append(field_alt_list)

    def get_schedule_field(self, target_field, repeat=1):

        _length   = self.check_length.to(u.hour).value
        _step     = self.check_step.to(u.hour).value
        _location = self.obs_location
        _alt_list    = self.obs_alt
        _target_list = target_field
        
        # get block time according to ra range
        block_time = (_target_list.ra_max - _target_list.ra_min) / (15.*u.deg) * u.hour
        self.block_time = block_time
        block_time = self.block_time

        field_start_time_list = []
        field_start_az_list   = []
        field_start_alt_list  = []
        nday = len(self.obs_date)
        dec_min = _target_list.dec_min
        dec_max = _target_list.dec_max
        dec_list = np.linspace(dec_min, dec_max, 2*nday/repeat, endpoint=False)
        ddec = dec_list[1] - dec_list[0]
        dec_list += ddec * 0.5
        dec_list = np.repeat(dec_list[None, :], repeat, axis=0).flatten()
        for oo,  obs_date in enumerate(self.obs_date):

            _alt = _alt_list[oo]

            _target0 = SkyCoord(ra = _target_list.ra, dec = dec_list[oo])
            _target1 = SkyCoord(ra = _target_list.ra, dec = dec_list[oo+nday])
            
            _obs_time = obs_date + np.arange(0, _length, _step) * u.hour
            rt = _get_risingtime_at_alt(_obs_time, _location, _alt, _target0)[0]
            st = _get_risingtime_at_alt(_obs_time, _location, _alt, _target1)[1]
            
            _obs_time = rt + (np.linspace(-1, 1, 100) * _step) * u.hour
            rt = _get_risingtime_at_alt(_obs_time, _location, _alt, _target0)[0]
            _raltaz = _target0.transform_to(AltAz(obstime=rt, location=_location))
            ralt, raz = _raltaz.alt.deg, _raltaz.az.deg
            # make sure the az range between -185 to 275
            if raz < -185.: raz += 360.
            if raz > 275. : raz -= 360.
            rt = rt - block_time[oo][0].to(u.hour) * 0.5

            _obs_time = st + (np.linspace(-1, 1, 100) * _step) * u.hour
            st = _get_risingtime_at_alt(_obs_time, _location, _alt, _target1)[1]
            _saltaz = _target1.transform_to(AltAz(obstime=st, location=_location))
            salt, saz = _saltaz.alt.deg, _saltaz.az.deg
            # make sure the az range between -185 to 275
            if saz < -185.: saz += 360.
            if saz > 275. : saz -= 360.
            st = st - block_time[oo][1].to(u.hour) * 0.5

            field_start_time_list.append([rt, st])
            field_start_az_list.append(  [raz  * u.deg, saz  * u.deg])
            field_start_alt_list.append( [ralt * u.deg, salt * u.deg])

        return field_start_time_list, field_start_az_list, field_start_alt_list

def make_obs_oneday(obs_az_range, obs_speed, obs_date_list, field, 
        alt_list, az_list, t_list, cal_list, obs_tot, location):

    status = ['Rising Block', 'Setting Block']

    _obs_info = []
    _msg = ''

    for oo, obs_date in enumerate(obs_date_list):

        for ii in range(2):

            _obs_az_range = obs_az_range[oo][ii]
            _obs_speed = obs_speed[oo][ii].to(u.arcmin/u.second).value
            _slewtime = (_obs_az_range/obs_speed[oo][ii]).to(u.second)
            _obs_tot  = obs_tot[oo][ii].to(u.second).value
            _slewnumb = int(round(_obs_tot / _slewtime.value / 2.))

            _msg += "# -- %14s --\n"%status[ii]
            _msg += "# ScanSpeed  %10.5f arcmin / s\n"%(_obs_speed)
            _msg += "# AZRange    %10.5f deg\n"%(_obs_az_range.to(u.deg).value)
            _msg += "# ScanTime   %10.5f s\n"%(_slewtime.value * 2)
            _msg += '# ScanNumb   %10d \n'%(_slewnumb)
            _msg += "# Calibrator %s + %s\n"%(cal_list[oo][ii * 2].name,
                                          cal_list[oo][ii * 2 + 1].name)

            _alt_0 = alt_list[oo][ii][0].to(u.deg).value
            _az_0 = az_list[oo][ii][0].to(u.deg).value
            if _az_0 > 275: _az_0 -= 360
            elif _az_0 < -185: _az_0 += 360
    
            syst_setup = 1.5 * u.min # for the system setup
            moving_speed = 2. * u.deg / u.second
            cal = cal_list[oo][ii * 2]
            cal_coord   = SkyCoord(cal.ra, cal.dec, frame='icrs')
            cal_altaz   = cal_coord.transform_to(AltAz(location=location,
                obstime=t_list[oo][ii][0] - cal.cal_time * u.min))
            #field_coord = SkyCoord(ra  = field.ra_min * u.deg, dec = field.dec)
            #sep_wgzmin = cal_coord.separation(field_coord)
            cal_az = cal_altaz.az.deg
            #if cal_az > 275: cal_az -= 360.
            #if cal_az < -185: cal_az += 360.
            if cal_az > 180: cal_az -= 360.
            if cal_az < -180: cal_az += 360.
            _msg += "# Init. Posi. Alt = %3.2f deg, Az = %3.2f deg\n"%(
                    cal_altaz.alt.deg, cal_az)
            sep_wgzmin = abs(cal_az - _az_0)
            cal_moving = (sep_wgzmin * u.deg / moving_speed).to(u.min)
            cal_moving += 10 * u.second
            cal_time   = cal.cal_time * u.min

            ahead_time = (syst_setup + cal_moving + cal_time).to(u.min)
            start_time = t_list[oo][ii][0]
            _start_time = t_list[oo][ii][0] - ahead_time
            _msg += "\n"
            _msg += "\tSuggest start time (%3.1f min ahead)\n"%( ahead_time.value)
            _msg += "\t%s [LST: %10.8f]\n"%(_start_time, 
                    _start_time.sidereal_time('apparent').value)
            _msg += "\t   system setup time %f min \n"%syst_setup.to(u.min).value +\
                    "\t  +calibration time  %f min \n"%cal_time.to(u.min).value +\
                    "\t  +moving from cal   %f min \n"%cal_moving.to(u.min).value
            _msg += "\n"

            _msg += "%s\t%s [LST: %10.8f] start\n"%('x', t_list[oo][ii][0],
                    t_list[oo][ii][0].sidereal_time('apparent').value)
            _msg += "\tAZ Alt: %10.8f, %10.8f\n"%(_az_0, _alt_0)
            _msg += "\tt_tot: %8.2f s [%8.2f min]\n"%(obs_tot[oo][ii].to(u.second).value,
                    obs_tot[oo][ii].to(u.min).value)
            _msg += "\t%s [LST: %10.8f] finish\n"%(t_list[oo][ii][-1],
                    t_list[oo][ii][-1].sidereal_time('apparent').value)

            cal = cal_list[oo][ii * 2 + 1]
            cal_coord   = SkyCoord(cal.ra, cal.dec, frame='icrs')
            cal_altaz   = cal_coord.transform_to(AltAz(location=location,
                obstime=t_list[oo][ii][-1]))
            #field_coord = SkyCoord(ra  = field.ra_max * u.deg, dec = field.dec)
            #sep_wgzmax = cal_coord.separation(field_coord)
            cal_az = cal_altaz.az.deg
            if cal_az > 275: cal_az -= 360.
            if cal_az < -185: cal_az += 360.
            sep_wgzmax = abs(cal_az - _az_0 - _obs_az_range.to(u.deg).value)
            cal_moving = (sep_wgzmax * u.deg / moving_speed).to(u.min)
            cal_moving += 10 * u.second
            cal_time   = cal.cal_time * u.min
            ahead_time = (cal_moving + cal_time).to(u.min)
            _msg += "\n"
            _msg += "\tFinishing time (%3.1f min added)\n"%ahead_time.value
            _msg += "\t%s [LST: %10.8f]\n"%(t_list[oo][ii][-1] + ahead_time, 
                    (t_list[oo][ii][-1] + ahead_time).sidereal_time('apparent').value)
            _msg += "\t  +calibration time  %f min \n"%cal_time.to(u.min).value +\
                    "\t  +moving to cal     %f min \n"%cal_moving.to(u.min).value
            _msg += '-' * 60
            _msg += '\n\n'

            _obs_info.append((
                status[ii], 
                '%s + %s'%( cal_list[oo][ii * 2].name, 
                                     cal_list[oo][ii * 2 + 1].name), 
                _obs_speed,
                _obs_tot,
                _obs_az_range.to(u.deg).value, 
                _slewtime.value, _slewnumb,
                '%s'%start_time, 
                start_time.sidereal_time('apparent').value,
                '%f %f'%(_az_0, _alt_0)))
    _obs_info = np.array(_obs_info, dtype = [
        ('Obs.', 'S20'), ('Calibrator', 'S15'), 
        ('Scan Speed\n [arcmin/s]', 'f4'),
        ('Total Time\n [s]', 'f4'),
        ('AZ Range\n [deg]', 'f4'), 
        ('Scan Time\n (One-way) [s]', 'f4'), 
        ('Scan No.\n (Two-way) [#]', 'i4'), 
        ('Start Time\n [UTC]', 'S25'), 
        ('Start Time\n [LST]', 'f4'), 
        ('Pointing\n (AZ Alt)', 'S20')])

    return _obs_info, _msg
    

def plot_sch_oneday(obstime_list, location, field_list, ax, obs_alt=None, cal_list=None):

    if obs_alt is not None:
        obs_alt = [x.to(u.deg).value for x in obs_alt]
    
    altaz_frame = AltAz(obstime=obstime_list, location=location)

    obs = Observer(location=location)
    print 'Sun rise time: %s'%(Time(obs.sun_rise_time(obstime_list[-1])).fits)
    print 'Sun set time: %s'%(Time(obs.sun_set_time(obstime_list[0])).fits)

    day = obstime_list[0]
    day.out_subfmt = 'date'
    #print day, day.sidereal_time('apparent')
    xx = mdates.date2num([t.datetime for t in obstime_list])

    alt_list = []

    sun_altaz = get_sun(obstime_list).transform_to(altaz_frame)
    yy = sun_altaz.alt.deg
    ax.plot(xx, yy, 'r-', label=r'${\rm SUN}$')

    moon_altaz = get_moon(obstime_list).transform_to(altaz_frame)
    yy = moon_altaz.alt.deg
    ax.plot(xx, yy, 'y-', label=r'${\rm MOON}$')

    if cal_list is not None:
        for cc, cal in enumerate(cal_list):
            cal_coord = SkyCoord(cal.ra, cal.dec, frame='icrs')
            cal_altaz = cal_coord.transform_to(altaz_frame)
            yy = cal_altaz.alt.deg
            ax.plot(xx, yy, color=_c_list[cc + len(field_list)], 
                    linestyle='--', label=r'${\rm %s}$'%cal.name)

            cal_altaz = cal_coord.transform_to(altaz_frame)
            sun_sep  = cal_altaz.separation(sun_altaz).to(u.deg).value
            moon_sep = cal_altaz.separation(moon_altaz).to(u.deg).value

            print 'Cal. %s moon sep: %5.2f-%5.2f deg; sun sep %5.2f-%5.2f deg.'%(
                    cal.name, moon_sep.min(), moon_sep.max(), 
                    sun_sep.min(), sun_sep.max())
    print '-'*10


    for ff, field in enumerate(field_list):
        #ra = np.mean(wgz_field[key][:2]) * u.deg + fcenter_shift[0] * u.deg
        #dec = np.mean(wgz_field[key][2:]) * u.deg + fcenter_shift[1] * u.deg
        ra  = field.ra
        dec = field.dec
        _wgz_altaz = SkyCoord(ra=ra, dec=dec).transform_to(altaz_frame)
        yy = _wgz_altaz.alt.deg
        #ax.plot(xx, yy, color=_c_list[ff], label=r'${\rm %s}$'%field.name)
        ax.plot(xx, yy, color=_c_list[ff], label=r'%s'%field.name)
    
        _wgz_altaz_min = SkyCoord(ra=field.ra_min, 
                                  dec=dec).transform_to(altaz_frame)
        _wgz_altaz_max = SkyCoord(ra=field.ra_max, 
                                  dec=dec).transform_to(altaz_frame)
        ax.fill_between(xx, y1 = _wgz_altaz_max.alt.deg, y2 = _wgz_altaz_min.alt.deg,
                        color = _c_list[ff], alpha=0.2)
        
        sun_sep = _wgz_altaz.separation(sun_altaz)
        #ax.plot(xx, sun_sep.deg, 'r--')
        moon_sep = _wgz_altaz.separation(moon_altaz)
        moon_sep_min = _wgz_altaz_min.separation(moon_altaz)
        moon_sep_max = _wgz_altaz_max.separation(moon_altaz)
        #ax.plot(xx, moon_sep.deg, 'y--')
    
    
        # find the time that targets are close to obs alt. 

        if obs_alt is not None:
            close_points = (_wgz_altaz.alt.deg >= obs_alt[0]).astype('int')
            close_points[1:] = close_points[1:] - close_points[:-1]
            close_points[0]  = 0
            r_points = close_points==1

            close_points = (_wgz_altaz.alt.deg >= obs_alt[1]).astype('int')
            close_points[1:] = close_points[1:] - close_points[:-1]
            close_points[0]  = 0
            s_points = close_points==-1

            ax.vlines(xx[r_points], -90, 90, colors='k', linestyles=':')
            ax.vlines(xx[s_points], -90, 90, colors='k', linestyles='--')

            ax.text(xx[r_points], 70, r'$%4.1f\,$'%sun_sep[r_points].deg, 
                    horizontalalignment='right', color='r')
            ax.text(xx[r_points], 60, r'$%4.1f\,[%4.1f, %4.1f]\,$'%(
                moon_sep[r_points].deg, 
                moon_sep_min[r_points].deg, 
                moon_sep_max[r_points].deg,), 
                horizontalalignment='right', color='y')

            ax.text(xx[s_points], 70, r'$%4.1f\,$'%sun_sep[s_points].deg, 
                    horizontalalignment='right', color='r')
            ax.text(xx[s_points], 60, r'$%4.1f\,[%4.1f, %4.1f]\,$'%(
                moon_sep[s_points].deg,
                moon_sep_min[s_points].deg,
                moon_sep_max[s_points].deg,),
                horizontalalignment='right', color='y')
    
    if obs_alt is not None:
        ax.hlines(obs_alt[0], xx[0], xx[-1], colors='0.5', linestyles='-', linewidth=2)
        ax.hlines(obs_alt[1], xx[0], xx[-1], colors='0.5', linestyles='-', linewidth=2)

    #date_format = mdates.DateFormatter('%m/%d %H:%M')
    date_format = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    #ax.xaxis.set_tick_params(rotation=90)
    ax.xaxis.set_tick_params(which='both', direction='in')
    #ax.legend(frameon=False, ncol=4)
    ax.set_ylabel('%s'%day)
    ax.set_ylim(ymin=0, ymax=90)
    ax.set_xlim(xmin=xx[0], xmax=xx[-1])

def _get_risingtime_at_alt(obstime, location, alt, target, return_idx=False):

    '''
    obstime : list of Time
    alt     : [rising_alt, setting_alt] observing alt
    target  : [ra, dec] of the target, in deg 
    '''

    idx = range(len(obstime))
    altaz_frame = AltAz(obstime=obstime, location=location)
    if return_idx: obstime=idx

    #ra, dec = target
    #ra  =  ra * u.deg
    #dec = dec * u.deg

    #target_altaz = SkyCoord(ra, dec).transform_to(altaz_frame)
    target_altaz = target.transform_to(altaz_frame)
    close_points = target_altaz.alt.deg > alt[0].to(u.deg).value
    close_points = close_points.astype('int')
    close_points[1:] = close_points[1:] - close_points[:-1]
    close_points[0]  = 0
    rising_point = (close_points == 1)
    if np.any(rising_point):
        rt = obstime[rising_point][0]
    else:
        rt = None

    close_points = target_altaz.alt.deg > alt[1].to(u.deg).value
    close_points = close_points.astype('int')
    close_points[1:] = close_points[1:] - close_points[:-1]
    close_points[0]  = 0
    setting_point = (close_points == -1)
    if np.any(setting_point):
        st = obstime[setting_point][0]
    else:
        st = None

    return rt, st



def make_obs(obs_az_range, obs_speed, obs_int, obs_schedul, field,
        location=None, az_inverse=False, cal=None):
    
    print '='*70

    print
    print "# Log.Lat.   %12.7f %12.7f" % (location.lon.deg, location.lat.deg)
    #print "# AZRange    %10.5f deg"%(obs_az_range.to(u.deg).value)
    
    print "# Field  %10s"%field.name
    print "# AZRange    %10.5f deg"%(obs_az_range.to(u.deg).value)
    print "# ScanSpeed  %10.5f arcmin / s"%(obs_speed.to(u.arcmin/u.second).value)
    print "# Int.Time   %10.5f s"%(obs_int.to(u.second).value)
    print "# SlewTime   %10.5f s"%((obs_az_range/obs_speed).to(u.second).value)
    #print '-'*70
    print

    for key in wgz_field.keys():
        #fig = plt.figure(figsize=(10, 4))
        #ax  = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        center_shift = {'R': 0.53, 'S': 0.47}

        for status in ['R', 'S']:

            ndays = len(obs_schedul[key][status])
            if ndays == 0: continue

            t_delay_pd = (obs_az_range * 2. / obs_speed).decompose() / float(ndays)
            #print t_delay_pd
            t_delay_pd = (np.arange(ndays) - 0.5 * (ndays - 1)) * t_delay_pd
            #t_delay_pd = t_delay_pd[::-1]

            dd = 0
            for o in obs_schedul[key][status]:
                _t, _alt, _az, obs_tot = o


                _t += t_delay_pd[dd]

                #print _t, _alt, _az

                _alt *= u.deg
                _az  *= u.deg

                dd += 1

                _az_space = obs_speed * obs_int
                _one_way_npoints = (obs_az_range / obs_speed / obs_int).decompose()
                _az_list  = np.arange(_one_way_npoints)\
                          - center_shift[status] * _one_way_npoints
                if az_inverse:
                    _az_list = _az_list[::-1]
                _az_list  = np.append(_az_list, _az_list[::-1])
                _az_list *= _az_space
                _az_list += _az

                obs_len = (obs_tot / obs_int).decompose()
                #_t_list = np.arange(obs_len) - 0.5 * obs_len
                _t_list = range(obs_len) - 0.5 * obs_len
                _t_list *= obs_int
                _t_list += _t

                _az_list = [_az_list[i%int(2*_one_way_npoints)] for i in range(obs_len)]

                _az_0 = _az_list[0].to(u.deg).value
                if _az_0 > 275: _az_0 -= 360
                elif _az_0 < -185: _az_0 += 360

                #if cal is not None:
                #    print "%s\t%si - %s"%(_t_list[0] - 12. * u.min, )
                ahead_time = 0 * u.min
                syst_setup = 1.5 * u.min # for the system setup
                moving_speed = 2. * u.deg / u.second
                cal_moving = 0 * u.min
                cal_time   = 0 * u.min
                #cal_moving = obs_schedul[key][status+'_cal_moving'] * u.min
                #if cal is not None:
                #    ahead_time += 14 * u.min + 1. * u.min
                #print "------\tSuggest start time (%3.1f min ahead)"%(
                #        ahead_time.value + cal_moving.value)
                #print "\t%s [LST: %10.8f]"%(_t_list[0] - ahead_time - cal_moving, 
                #(_t_list[0] - ahead_time - cal_moving).sidereal_time('apparent').value)
                if cal is not None:
                    cal_key = obs_schedul[key][status+"_cal"][0]
                    cal_coord = SkyCoord(cal[cal_key][0], cal[cal_key][1],
                                         frame='icrs')
                    #altaz_frame = AltAz(obstime=_t_list[0] - ahead_time,
                    #                    location=location)
                    altaz_frame = AltAz(obstime=_t_list[0], location=location)
                    _wgz_altaz_min = SkyCoord(ra=wgz_field[key][0] * u.deg,
                            dec=wgz_field[key][3] * u.deg).transform_to(altaz_frame)
                    moon_altaz = get_moon(_t_list[0]).transform_to(altaz_frame)
                    moon_radec = get_moon(_t_list[0])
                    cal_altaz  = cal_coord.transform_to(altaz_frame)
                    cal_az  = cal_altaz.az.deg
                    if cal_az > 275: cal_az -= 360.
                    elif cal_az < -185: cal_az += 360.
                    cal_alt = cal_altaz.alt.deg
                    sep_moon   = cal_altaz.separation(moon_altaz)
                    sep_wgzmin = cal_altaz.separation(_wgz_altaz_min)
                    cal_moving = (sep_wgzmin.deg * u.deg / moving_speed).to(u.min)
                    cal_moving += 10 * u.second
                    cal_time   = obs_schedul[key][status+'_cal_time'] * u.min

                ahead_time = (syst_setup + cal_moving + cal_time).to(u.min)
                print '-' * 60
                print "\tSuggest start time (%3.1f min ahead)"%(
                        ahead_time.value)
                print "\t%s [LST: %10.8f]"%(_t_list[0] - ahead_time, 
                (_t_list[0] - ahead_time).sidereal_time('apparent').value)
                print "\t   system setup time %f min \n"%syst_setup.to(u.min).value +\
                      "\t  +calibration time  %f min \n"%cal_time.to(u.min).value +\
                      "\t  +moving from cal   %f min \n"%cal_moving.to(u.min).value
                    
                if cal is not None:
                    print
                    print "C\t%s [RA Dec: %10.8f, %10.8f]"%(cal_key,
                                cal_coord.ra.deg, cal_coord.dec.deg)
                    print " \t%s [AZ Alt: %10.8f, %10.8f]"%(cal_key, cal_az, cal_alt)
                    print "\tStart point Sep. %10.8f "%(sep_wgzmin.deg)
                    print "\tMoon Sep. %10.8f [RA: %10.8f Dec: %10.8f]"%(sep_moon.deg,
                            moon_radec.ra.deg, moon_radec.dec.deg)
                print
                print "%s\t%s [LST: %10.8f] start"%(status, _t_list[0],
                                    _t_list[0].sidereal_time('apparent').value)
                print "\tAZ Alt: %10.8f, %10.8f"%(_az_0, _alt.value)
                print "\tt_tot: %8.2f s [%8.2f min]"%(obs_tot.to(u.second).value,
                        obs_tot.to(u.min).value)
                print "\t%s [LST: %10.8f] finish"%(_t_list[-1],
                                    _t_list[-1].sidereal_time('apparent').value)
                if cal is not None:
                    altaz_frame = AltAz(obstime=_t_list[-1], location=location)
                    cal_key = obs_schedul[key][status+"_cal"][1]
                    cal_coord = SkyCoord(cal[cal_key][0], cal[cal_key][1],
                                         frame='icrs')
                    _wgz_altaz_max = SkyCoord(ra=wgz_field[key][1] * u.deg,
                            dec=wgz_field[key][2] * u.deg).transform_to(altaz_frame)
                    cal_altaz  = cal_coord.transform_to(altaz_frame)
                    sep_wgzmax = cal_altaz.separation(_wgz_altaz_max)
                    cal_az  = cal_altaz.az.deg
                    if cal_az > 275: cal_az -= 360.
                    elif cal_az < -185: cal_az += 360.
                    cal_alt = cal_altaz.alt.deg
                    cal_moving = (sep_wgzmax.deg * u.deg / moving_speed).to(u.min)
                    cal_moving += 10 * u.second
                    cal_time   = obs_schedul[key][status+'_cal_time'] * u.min

                if cal is not None:
                    print
                    print "C\t%s [RA Dec: %10.8f, %10.8f]"%(cal_key,
                                cal_coord.ra.deg, cal_coord.dec.deg)
                    print " \t%s [AZ Alt: %10.8f, %10.8f]"%(cal_key, cal_az, cal_alt)
                    print "\tFinish point Sep. %10.8f "%(sep_wgzmax.deg)

                ahead_time = (cal_moving + cal_time).to(u.min)
                print
                print "\tFinishing time (%3.1f min added)"%ahead_time.value
                print "\t%s [LST: %10.8f]"%(_t_list[-1] + ahead_time, 
                        (_t_list[-1] + ahead_time).sidereal_time('apparent').value)
                print "\t   system setup time %f min \n"%syst_setup.to(u.min).value +\
                      "\t  +calibration time  %f min \n"%cal_time.to(u.min).value +\
                      "\t  +moving to cal     %f min \n"%cal_moving.to(u.min).value
                print '-' * 60
                print

                #pp = SkyCoord(alt=_alt, az=_az_list, frame='altaz',
                #         location=location, obstime=_t_list)
                #pp = pp.transform_to('icrs')

                #radec_list = np.array([[_pp.ra.deg, _pp.dec.deg] for _pp in pp])

                #if key in ['00h', '01h']:
                #    radec_list[:,0][radec_list[:,0]>180] -= 360

def get_az_raster_range(patch, alt_list, obs_location):
    dec_min = patch.dec_min * np.pi/180.
    dec_max = patch.dec_max * np.pi/180.

    lat = obs_location.lat.rad

    alt_list = np.array(alt_list)
    alt_list = alt_list * np.pi/180.

    az_1 = np.arccos( (np.sin(dec_max) - np.sin(alt_list) * np.sin(lat)) / (np.cos(alt_list) * np.cos(lat)) )
    az_2 = np.arccos( (np.sin(dec_min) - np.sin(alt_list) * np.sin(lat)) / (np.cos(alt_list) * np.cos(lat)) )

    az_diff = np.abs(az_1 - az_2)

    for ii, _x in enumerate(az_diff):
        if not np.isfinite(_x):
            print "Warning: Alt = %3.1f is too high to cover the full Dec range"%\
            (alt_list[ii] * 180./np.pi)

    az_diff *= 180./np.pi
    return az_diff

def get_scanning_speed(alt_list, speed = 5.):

    return speed / np.cos(np.array(alt_list) * np.pi / 180.)
    
if __name__ == "__main__":

    pass
