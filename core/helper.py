import datetime
import logging
import math
import os
import re
import sys

from datetime import timedelta

import georinex as gr
import numpy as np
import scipy.interpolate as interp
import dateparser

import settings as settings

from app.downloads import DownloadDCB, DownloadOrbit, DownloadGlonassChannel
from app.parser import ParserChannels, ParserRinexChannels, ParserOrbit, ParserDCB


class Utils:
    """
    Python utils, converting and dealing with variable types
    """

    @staticmethod
    def array_timestamp_to_datetime(array_timestamp):
        """
        Convert an array with string timestamp to an array of datetimes

        :param array_timestamp: Array with string date in timestamp format
        :return: Array with datetime format elements
        """
        array = []

        for item in array_timestamp.values:
            array.append(dateparser.parse(str(item)))

        return array

    @staticmethod
    def array_dict_to_array(array_dict, key):
        """
        Get all values from a specific key in a dictionary array, and transform in a simple array

        :param array_dict: Array of dictionaries
        :param key: Interest key, based on this key that the values in array_dict will be extracted
        :return: Array, with values based on the key passed on
        """
        array = []

        for item in array_dict:
            array.append(item[key])

        return array

    @staticmethod
    def rinex_str_date_to_datetime(string_date):
        """
        From rinex string timestamps, return the right datetime
        :param string_date:
        :return:
        """
        m = re.search(settings.REGEX_RINEX_DATE, string_date)
        year = int(m.group(1))
        month = int(m.group(2))
        day = int(m.group(3))
        hour = int(m.group(4))
        minute = int(m.group(5))
        second = int(float(m.group(6)))

        return datetime.datetime(year, month, day, hour, minute, second)

    @staticmethod
    def interpolate_orbit(orbit, hdr):
        """
        Interpolate X, Y, Z values as sampling on the Rinex file. Thus, both measures have the same precision

        :param orbit: The python orbit object less precise than rinex
        :param hdr: Rinex header, to extract rinex precision, first and last acquisition
        :return: The python orbit object with the same precise than rinex, with interpolated values
        """
        orbit_interpolated = {}
        rinex_last_obs_time = Utils.rinex_str_date_to_datetime(hdr['TIME OF LAST OBS'])

        for prn, data in orbit.items():
            date_array_in_orbit = np.array([])
            date_array_interpolated = np.array([])
            x_aux = np.array([])
            y_aux = np.array([])
            z_aux = np.array([])

            for item in data:
                if item['date'] > rinex_last_obs_time:
                    continue

                date_array_in_orbit = np.append(date_array_in_orbit, item['date'])
                x_aux = np.append(x_aux, item['x'])
                y_aux = np.append(y_aux, item['y'])
                z_aux = np.append(z_aux, item['z'])

            date_array_in_orbit = np.append(date_array_in_orbit, rinex_last_obs_time)

            aux_date = date_array_in_orbit[0]
            date_array_interpolated = np.append(date_array_interpolated, aux_date)
            while date_array_in_orbit[-1] != aux_date:
                aux_date += timedelta(seconds=float(hdr['interval']))
                date_array_interpolated = np.append(date_array_interpolated, aux_date)

            x_interp = interp.CubicSpline(np.arange(x_aux.size), x_aux)
            x_interp = x_interp(np.linspace(0, x_aux.size - 1, date_array_interpolated.size))

            y_interp = interp.CubicSpline(np.arange(y_aux.size), y_aux)
            y_interp = y_interp(np.linspace(0, y_aux.size - 1, date_array_interpolated.size))

            z_interp = interp.CubicSpline(np.arange(z_aux.size), z_aux)
            z_interp = z_interp(np.linspace(0, z_aux.size - 1, date_array_interpolated.size))

            orbit_interpolated[prn] = []
            for i in range(len(date_array_interpolated)):
                orbit_interpolated[prn].append({'date': date_array_interpolated[i], 'x': x_interp[i],
                                                'y': y_interp[i], 'z': z_interp[i]})

        return orbit_interpolated

    @staticmethod
    def convert_dcb_ns_to_meter(dcb):
        """
        DCB files originally store values in nanoseconds unit. In order to use those values to remove ou correct another
        files or calculation. The units has to be the same. Here, these values are converted to meters unit

        :param dcb: The original DCB parsed file, with nanoseconds values
        :return: The updated DCB parsed file, with the values of bias updated to meters
        """
        dcb_ns = {}

        if not dcb:
            logging.info(">>>> DCB is empty. Skipping convertion!")
            return

        for key, dcb_type in dcb.items():
            dcb_ns[key] = {}
            dcb_aux = {}
            for prn, value in dcb_type.items():
                ns_to_meter = settings.C * pow(10, -9)
                new_value = [x * ns_to_meter for x in value]
                dcb_aux[prn] = [new_value[0], value[1]]

            dcb_ns[key].update(dcb_aux)

        return dcb_ns

    @staticmethod
    def convert_dcb_ns_to_tecu(dcb, factor_glonass):
        """
        DCB files originally store values in nanoseconds unit. In order to use those values to remove ou correct another
        files or calculation. The units has to be the same. Here, these values are converted to TEC unit (TECU)

        :param dcb: The original DCB parsed file, with nanoseconds values
        :param factor_glonass: Python object with GLONASS factors
        :return: The updated DCB parsed file, with the values of bias updated to TEC Unit
        """
        dcb_tecu = {}

        for key, dcb_type in dcb.items():
            dcb_tecu[key] = {}
            dcb_aux = {}
            for prn, value in dcb_type.items():
                const_str = prn[0:1]
                prn_str = prn[1:]

                ns_to_meter = settings.C * pow(10, -9)

                if const_str == 'G':
                    # TODO: conferir com Haroldo/Cosme sobre o uso do factor para converter em TECU
                    # factor_gps = (f1 - f2) / (f1 + f2) / settings.C
                    f1 = settings.FREQUENCIES['G']['1']
                    f2 = settings.FREQUENCIES['G']['2']
                    factor_gps = pow(f1 * f2, 2) / (f1 + f2) / (f1 - f2) / settings.A / settings.TECU
                    new_value = [(x * ns_to_meter * factor_gps) for x in value]
                    dcb_aux[prn] = [new_value[0], value[1]]
                elif const_str == 'R' and prn_str in factor_glonass:
                    new_value = [(x * ns_to_meter * factor_glonass[prn_str][2]) for x in value]
                    dcb_aux[prn] = [new_value[0], value[1]]

            dcb_tecu[key].update(dcb_aux)

        return dcb_tecu

    @staticmethod
    def check_availability(factor_glonass, dcb, constellation, prn):
        """
        Sometimes, the GLONASS factor or DCB array values, might not be available for some specific PRN. Thus, this
        method apply default values for both variables for such cases

        :param factor_glonass: GLONASS factor, with factors values for GLONASS constellation
        :param dcb: Dict object with bias of C1-P1 and P2-P1
        :param constellation: The constellation verified: 'G', 'R', etc
        :param prn: The satellite code (PRN) verified: '01', '02', etc
        :return: In case the data are not available, it returns the default values. For instance, the value of 0 for DCB
        indicates it has no compensation under a respective equation
        """
        dcb_compensate = 0
        f1 = settings.FREQUENCIES[constellation]['1']
        f2 = settings.FREQUENCIES[constellation]['2']

        if constellation == 'G':
            factor = pow(f1 * f2, 2) / (f1 + f2) / (f1 - f2) / settings.A / settings.TECU
        elif constellation == 'R':
            if prn not in factor_glonass:
                logging.info(">>>>>> Relative TEC for {} constellation, was not compensate due to "
                             "the lack of factor GLONASS for PRN {}".format(constellation, prn))
            else:
                factor = factor_glonass[prn][2]

        key_in_dcb = constellation + prn

        if key_in_dcb not in dcb['P1-P2']:
            logging.info(">>>>>> Relative TEC for {} constellation, was not compensate due to "
                         "the lack of DCB for PRN {}".format(constellation, prn))
        else:
            dcb_compensate = dcb['P1-P2'][key_in_dcb][0]

        return factor, dcb_compensate


class InputFiles:
    """
    Complementary methods to download all the input files needed for TEC and Bias estimation procedures
    """

    def _setup_rinex_name(self, rinex_folder, rinex_name):
        """
        Test if rinex name is in the old fashion way or in another formats. In case the format is newer or older, the
        method will always return the values needed
            Example of formats name accept:
                2.11: ALMA1520.18O
                3.03: ALMA00BRA_R_20181520000_01D_30S_MO.rnx

        :param rinex_folder: Absolute path to the rinexes folder
        :param rinex_name: String rinex filename
        :return: Returns the String rinex absolute path, and the year, month and doy related
        """
        if rinex_name.endswith(".rnx"):
            day_i, day_f = 16, 19
            year_i, year_f, year_type = 12, 16, "%Y"
            extens = "[rR][nN][xX]$"
        elif bool(re.match("[oO]$", rinex_name[-1:])):
            day_i, day_f = 4, 7
            year_i, year_f, year_type = -3, -1, "%y"
            extens = "[\\d]{2}[oO]$"
        else:
            logging.error(">>>> Error during rinex file reading. Check it and try again!")
            sys.exit()

        if len(rinex_name) == 0:
            logging.error('>> Something wrong with parameter \'rinex_name\'!. Empty name!')
            sys.exit()
        elif not rinex_name[0:4].isalpha():
            logging.error('>> Something wrong with parameter \'rinex_name\'!. IAGA code not well format!')
            sys.exit()
        elif not rinex_name[day_i:day_f].isdigit() or int(rinex_name[day_i:day_f]) > 366:
            logging.error('>> Something wrong with parameter \'rinex_name\'!. Invalid day of the year!')
            sys.exit()
        elif not bool(re.match(extens, rinex_name[-3:])):
            logging.error('>> Something wrong with parameter \'rinex_name\'!. Wrong extension or not well format!')
            sys.exit()

        path = os.path.join(rinex_folder, rinex_name)
        doy = rinex_name[day_i:day_f]
        year = datetime.datetime.strptime(rinex_name[year_i:year_f], year_type).strftime('%Y')
        month = datetime.datetime.strptime(doy, '%j').strftime('%m')

        return path, year, month, doy

    def _setup_file_and_download(self, **kwargs):
        """
        Prepare the inputs, such as the correct name and paths before actually download it. The files, in this case,
        can only be DCB or Orbit

        :param year: String year [YYYY]
        :param month: String month [mm]
        :param day: String day [dd]
        :param rinex_hdr: Current rinex header. If there is at least one constellation without L1P column, the DCB is
        then downloaded. Otherwise, the L1P measures are already sufficient and the download is skipped
        :param file_type: The file type to be downloaded. It might be 'DCB' or 'Orbit'
        :return: Python object with a parsed file
        """
        obj = {}
        year = kwargs.get('year')
        month = kwargs.get('month')
        day = kwargs.get('day')
        rinex_hdr = kwargs.get('rinex_hdr')
        file_type = kwargs.get('file_type')

        requiried_version = str(settings.REQUIRED_VERSION)
        col_var = settings.COLUMNS_IN_RINEX[requiried_version]

        if file_type == "DCB":
            l1p_present = True

            for item in rinex_hdr['fields']:
                if item not in settings.CONSTELATIONS:
                    continue

                if col_var[item]['P1'] not in rinex_hdr['fields'][item]:
                    l1p_present = False
                    break

            if l1p_present:
                logging.info(">>>> P1 present in all columns!")

            dcb = DownloadDCB(year, day)
            dcb.download()

            dcb_parsed = ParserDCB(dcb.file_uncompressed)
            dcb_parsed.parser()

            obj = dcb_parsed.parsed

        elif file_type == "Orbit":

            logging.info(">> Searching " + file_type + " file...")

            orbit = DownloadOrbit(year, month, day)
            orbit.download()

            orbit_parsed = ParserOrbit(orbit.file_uncompressed)
            orbit_parsed.parser()

            return orbit_parsed.parsed

        return obj

    def _updating_c1_to_p1(self, hdr, obs, dcb_m):
        """
        Update the rinex C1 measures, converting it to P1 values, corresponding to P1 = C1 + DCB_P1-C1, and DCB in
        meters unit

        :param hdr: Current rinex header
        :param obs: The current rinex measures to be updated
        :param dcb_m: DCB values in meter unit
        :return: The updated current rinex with P1 values calculated through C1 and DCB
        """
        if not dcb_m:
            logging.info(">>>> DCB is empty. Skipping conversion!")
            return

        requiried_version = str(settings.REQUIRED_VERSION)
        cols_var = settings.COLUMNS_IN_RINEX[requiried_version]

        for prn in obs.sv.values:
            if cols_var[prn[0:1]]['P1'] not in hdr['fields'][prn[0:1]] and \
                    cols_var[prn[0:1]]['C1'] not in hdr['fields'][prn[0:1]]:
                logging.error(">>>> Columns {} and {} not in rinex for {}...".format(cols_var[prn[0:1]]['P1'],
                                                                                     cols_var[prn[0:1]]['C1'],
                                                                                     prn[0:1]))
                sys.exit()

            if cols_var[prn[0:1]]['P1'] not in hdr['fields'][prn[0:1]] and \
                    cols_var[prn[0:1]]['C1'] in hdr['fields'][prn[0:1]]:
                obs[cols_var[prn[0:1]]['C1']].sel(sv=prn).values -= dcb_m['C1-P1'][prn][0]

        return obs

    def _which_cols_to_load(self):
        """
        The rinex file is an extensive file, sometimes, with a lot of measures that are not interesting for this present
        work. Said that, in this method is selected only the columns used during the EMBRACE TEC and Bias estimation

        :return: The columns to load in rinex
        """
        columns_to_be_load = []
        requiried_version = str(settings.REQUIRED_VERSION)
        dict_list_cols = list(settings.COLUMNS_IN_RINEX[requiried_version].values())

        for constellation in dict_list_cols:
            for item in constellation.values():
                if item not in columns_to_be_load:
                    columns_to_be_load.append(item)

        return columns_to_be_load

    def _prepare_rinex(self, complete_path, **kwargs):
        """
        Prepare the rinexs, such as checking if it is a valid file and if it has all the information need for the
        calculation

        :param complete_path: The absolute path to the rinex
        :param kwargs:
            is_partial: boolean that specify the moment rinex is read.
                Full-day: is_partial = False
                Partial: is_partial = True
            initial: ...
            final: ...
        :return: The Python objects after to successful read the prev and rinex, it includes the header
        and measures of both files
        """
        is_partial = kwargs.get('is_partial')

        if not os.path.isfile(complete_path):
            logging.error(">>>> Rinex " + complete_path + " does not exist. Process stopped!")
            sys.exit()

        columns_to_be_load = self._which_cols_to_load()

        logging.info(">> Reading rinex: " + complete_path)
        hdr = gr.rinexheader(complete_path)
        if hdr['version'] < settings.REQUIRED_VERSION:
            logging.error(">>>> Version {} required. Current version {}. "
                          "Process stopped!".format(settings.REQUIRED_VERSION, hdr['version']))
            sys.exit()

        if not is_partial:
            obs = gr.load(complete_path, meas=columns_to_be_load, use=settings.CONSTELATIONS)
        else:
            initial = kwargs.get('initial')
            final = kwargs.get('final')
            obs = gr.load(complete_path, meas=columns_to_be_load, use=settings.CONSTELATIONS, tlim=[initial, final])

        logging.info(">>>> Rinex version {}. Constellations: {}".format(hdr['version'], settings.CONSTELATIONS))

        return hdr, obs

    def _prepare_factor(self, hdr, year, day, month):
        """
        The factors is used as a weight during the first estimate (relative TEC). For the GLONASS constellation, the
        factors are selective, per PRN, then, this values can be download from the own rinex or a URL

        :param hdr: The current rinex header, which brings the information of GLONASS channels
        :param year: Year of the current rinex
        :param month: Month of the current rinex
        :param day: Day of the current rinex
        :return: The GLONASS factors, already multiply by the frequencies F1 and F2
            GLONASS factor Python object. Format example:
                    factor_glonass = {
                                '01': VALUE_1, VALUE_2, VALUE_3
                                '02': VALUE_1, VALUE_2, VALUE_3
                                '03': VALUE_1, VALUE_2, VALUE_3
                                '04': VALUE_1, VALUE_2, VALUE_3
                                '05': ...
                        }
        """
        logging.info(">> Downloading GLONASS channels in " + settings.URL_GLONASS_CHANNELS
                     + " for pseudorange calculation...")

        if "GLONASS SLOT / FRQ #" not in hdr:
            glonnas_channel = DownloadGlonassChannel(year, day, month)
            glonnas_channel.download()

            glonass_channels_parser = ParserChannels(glonnas_channel.file_uncompressed)
            glonass_channels_parser.parser()
        else:
            glonass_channels_parser = ParserRinexChannels(hdr['GLONASS SLOT / FRQ #'])
            glonass_channels_parser.parser()

        return glonass_channels_parser.parsed

    def _prepare_orbit(self, hdr, year, month, day):
        """
        Construct the Orbit file and Python object. The file has to be downloaded from a remote repository. As soon it
        is guaranteed, the Python object can then be constructed.

        :param hdr: The current rinex header
        :param year: Year of the current rinex
        :param month: Month of the current rinex
        :param day: Day of the current rinex
        :return: The Python orbit object, with the right location of every satellite
                 Orbit Python object. Format example:
                    orbit = {
                        'G01': [{'date': VALUE, 'x': VALUE, 'y': VALUE, 'z': VALUE}, {...}}
                        'G02': [{'date': VALUE, 'x': VALUE, 'y': VALUE, 'z': VALUE}, {...}}
                        'G03': [{'date': VALUE, 'x': VALUE, 'y': VALUE, 'z': VALUE}, {...}}
                        ...
                        'R01': [{'date': VALUE, 'x': VALUE, 'y': VALUE, 'z': VALUE}, {...}}
                        'R02': [{'date': VALUE, 'x': VALUE, 'y': VALUE, 'z': VALUE}, {...}}
                        'R03': [{'date': VALUE, 'x': VALUE, 'y': VALUE, 'z': VALUE}, {...}}
                        ...
                    }
        """
        logging.info(">> Downloading Orbit files...")
        orbit = self._setup_file_and_download(year=year, month=month, day=day, rinex_hdr=hdr, file_type='Orbit')

        logging.info(">> Checking precision in time of both files, rinex and orbit...")
        if not math.isnan(hdr['interval']):
            rinex_interval = float(hdr['interval'])

        orbit_prn = list(orbit.keys())
        orbit_diff = orbit[orbit_prn[0]][1]['date'] - orbit[orbit_prn[0]][0]['date']
        orbit_interval = orbit_diff.total_seconds()

        if rinex_interval < orbit_interval:
            logging.info(">>>> Orbit is less precise ({}) than rinex ({}) in seconds. Interpolation will be applied!".
                         format(orbit_interval, rinex_interval))

            orbit_aux = Utils.interpolate_orbit(orbit, hdr)
            orbit.update(orbit_aux)
        elif rinex_interval == orbit_interval:
            logging.info(">>>> Rinex and orbit with same time precision, {} and {} in seconds. "
                         "Interpolation is not required!".format(rinex_interval, orbit_interval))

        return orbit

    def _prepare_dcb(self, hdr, year, month, day):
        """
        Construct the DCB file and Python object. The file has to be downloaded from a remote repository. As soon it
        is guaranteed, the Python object can then be constructed. Finally, the DCB values are converted from nanosecs to
        TEC units

        :param hdr: The current rinex header
        :param year: Year of the current rinex
        :param month: Month of the current rinex
        :param day: Day of the current rinex
        :return: The Python DCB object, with the converted TEC units satellites errors
            DCB Python object. Format example:
                    dcb = {
                        'C1-P1':
                            {'G01': [BIAS_VALUE, STD_VALUE],
                            'G02': [BIAS_VALUE, STD_VALUE],
                            'G03': [BIAS_VALUE, STD_VALUE],
                            ...
                            'G32': [BIAS_VALUE, STD_VALUE],
                            'R01': [BIAS_VALUE, STD_VALUE],
                            'R02': [BIAS_VALUE, STD_VALUE],
                            ...
                            'R24': [BIAS_VALUE, STD_VALUE]}
                        'P1-P2':
                            {'G01': [BIAS_VALUE, STD_VALUE],
                            'G02': [BIAS_VALUE, STD_VALUE],
                            'G03': [BIAS_VALUE, STD_VALUE],
                            ...
                            'G32': [BIAS_VALUE, STD_VALUE],
                            'R01': [BIAS_VALUE, STD_VALUE],
                            'R02': [BIAS_VALUE, STD_VALUE],
                            ...
                            'R24': [BIAS_VALUE, STD_VALUE]}
                    }
        """
        logging.info(">> Downloading DCB files...")
        dcb = self._setup_file_and_download(year=year, month=month, day=day, rinex_hdr=hdr, file_type='DCB')

        return dcb

    def prepare_rinex_partial(self, folder, file):
        """

        :param folder:
        :param file:
        :return:
        """
        path, year, month, doy = self._setup_rinex_name(folder, file)

        initial_date = datetime.datetime(int(year), 1, 1, settings.INITIAL_HOUR_RECALC_BIAS, 0) + \
                       datetime.timedelta(int(doy) - 1)
        final_date = datetime.datetime(int(year), 1, 1, settings.FINAL_HOUR_RECALC_BIAS, 0) + \
                     datetime.timedelta(int(doy) - 1)

        hdr, obs = self._prepare_rinex(path, is_partial=True, initial=initial_date, final=final_date)

        return hdr, obs

    def prepare_inputs(self, folder, file):
        """
        Prepare the inputs needed to the TEC and bias receptor estimation. It includes the current and previous day rinex,
        the orbit file, the DCBs, and the GLONASS channels

        :param folder: Absolute path to the rinex file
        :param file_prev: The previous rinex according to the current day
        :param file: Name of the file in order to extract the corresponding year and month values: ssssdddd.yyo
        :return: The rinex header and measures objects, orbit object, dcb object, and factor_glonass object
        """
        path, year, month, doy = self._setup_rinex_name(folder, file)

        hdr, obs = self._prepare_rinex(path, is_partial=False)
        factor_glonass = self._prepare_factor(hdr, year, month, doy)
        orbit = self._prepare_orbit(hdr, year, month, doy)
        dcb = self._prepare_dcb(hdr, year, month, doy)

        logging.info(">> Converting DCB nanoseconds in meter unit...")
        dcb_m = Utils.convert_dcb_ns_to_meter(dcb)

        logging.info(">> Converting C1 values to correspond P1 in rinex (P1 = C1 + DCB_P1-C1[meters])...")
        obs = self._updating_c1_to_p1(hdr, obs, dcb_m)

        return hdr, obs, orbit, dcb, factor_glonass


class Geodesy:
    """
    Complementary methods to calculate the right point of estimation, respect to receiver and satellites
    """

    def sub_ion_point(self, altitude, x_sat, y_sat, z_sat, x_rec, y_rec, z_rec):
        """
        Procedure to calculate the sub-ionospheric point through satellite and receiver positions

        :param altitude: Satellite altitude
        :param x_sat: X axis satellite position array
        :param y_sat: Y axis satellite position array
        :param z_sat: Z axis satellite position array
        :param x_rec: X axis receiver position value
        :param y_rec: Y axis receiver position value
        :param z_rec: Z axis receiver position value
        :return: Sub-ionospheric position array, corresponding to Lat, Long and Altitude, respectively
        """
        sub_ion_x = []
        sub_ion_y = []
        sub_ion_z = []

        for i, item in enumerate(x_sat):
            diff_x = x_sat[i] - x_rec
            diff_y = y_sat[i] - y_rec
            diff_z = z_sat[i] - z_rec

            phi = (x_sat[i]) * y_rec - (y_sat[i] * x_rec)
            theta = (x_sat[i] * z_rec) - (z_sat[i] * x_rec)

            sqr_sum = math.pow(diff_y, 2) + math.pow(diff_z, 2) + math.pow(diff_x, 2)
            factor = math.pow(diff_y * phi + diff_z * theta, 2) - sqr_sum * (math.pow(phi, 2) + math.pow(theta, 2) -
                                                                             math.pow(diff_x, 2) * math.pow(altitude,
                                                                                                            2))
            factor = math.sqrt(factor)

            dum = (-1 * (diff_y * phi + diff_z * theta) - factor) / sqr_sum

            if ((x_sat[i] - dum) * (x_rec - dum)) > 0:
                dum = (-1 * (diff_y * phi + diff_z * theta) + factor) / sqr_sum

            sub_ion_x.append(dum)
            sub_ion_y.append((diff_y * dum + phi) / diff_x)
            sub_ion_z.append((diff_z * dum + theta) / diff_x)

        return sub_ion_x, sub_ion_y, sub_ion_z

    def calc_zenital_angle(self, array_x_sat, array_y_sat, array_z_sat, x_rec, y_rec, z_rec):
        """
        Procedure to calculate the zenital angle through satellite and receiver positions

        :param array_x_sat: X axis satellite position array
        :param array_y_sat: Y axis satellite position array
        :param array_z_sat: Z axis satellite position array
        :param x_rec: X axis receiver position value
        :param y_rec: Y axis receiver position value
        :param z_rec: Z axis receiver position value
        :return: The zenital angle based on the satellite and receiver positions
        """
        x_sat = np.array(array_x_sat)
        y_sat = np.array(array_y_sat)
        z_sat = np.array(array_z_sat)

        ang = ((x_sat - x_rec) * x_rec) + ((y_sat - y_rec) * y_rec) + ((z_sat - z_rec) * z_rec)
        aux1 = pow(x_sat - x_rec, 2) + pow(y_sat - y_rec, 2) + pow(z_sat - z_rec, 2)
        aux2 = pow(x_rec, 2) + pow(y_rec, 2) + pow(z_rec, 2)
        aux1_sqrt = np.sqrt(aux1)
        aux2_sqrt = np.sqrt(aux2)

        result_ang_zen = ang / aux1_sqrt / aux2_sqrt
        ang_zen = np.arccos(result_ang_zen)
        ang_zen *= 180 / math.pi

        return ang_zen

    def car_2_pol(self, array_x, array_y, array_z):
        """
        Convert the normal cartesian projection in polar projection

        :param array_x: Subionospheric coordinate in x
        :param array_y: Subionospheric coordinate in y
        :param array_z: Subionospheric coordinate in z
        :return: The converted triplice coordination array
        """
        array_x_pol = []
        array_y_pol = []
        array_z_pol = []

        for i in range(len(array_x)):
            x = array_x[i]
            y = array_y[i]
            z = array_z[i]

            euclidian_distance = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
            e = (2 - (1 / settings.ELLIPTICITY)) / settings.ELLIPTICITY
            x_pol = 0.0
            y_pol = 0.0
            z_pol = 0.0

            if euclidian_distance != 0:
                x_pol = math.atan(z / euclidian_distance)

                for i in range(7):
                    n = settings.RADIUS_EQUATOR / math.sqrt(1 - e * math.pow(math.sin(z_pol), 2))
                    x_pol = math.atan(z / (euclidian_distance - (e * n * math.cos(x_pol))))

                if x != 0:
                    y_pol = math.atan(y / x)
                    if x < 0:
                        y_pol = y_pol + math.pi
                else:
                    y_pol = y / abs(y) * (math.pi / 2)

                z_pol = euclidian_distance / math.cos(x_pol) - n

            else:
                logging.info(">>>>>> Square root equal zero!")

            x_pol = x_pol * 180 / math.pi
            y_pol = y_pol * 180 / math.pi

            aux = y_pol
            if y_pol > 180.0:
                aux = y_pol - 180
                aux = -1 * (180 - aux)
            elif y_pol < -180.0:
                aux = y_pol + 180
                aux = 180 + aux

            y_pol = aux

            array_x_pol.append(x_pol)
            array_y_pol.append(y_pol)
            array_z_pol.append(z_pol)

        return array_x_pol, array_y_pol, array_z_pol

