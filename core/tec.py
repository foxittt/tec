import time
import os
import logging
import datetime
import csv

import core.settings as settings
import core.helper as helper
import core.estimate as calc
import core.cycleslip.cycleslip as cs


class TEC:
    """
    Class with the main calls for TEC estimation (by station)
    """
    def simplify_dict(self, dict_var, list_keys):
        """
        According to args, select and save in dict_var, specific keys among many other

        :param dict_var: Dict with all the TEC estimates
        :param list_keys: Brings the right keys to simplify the dict
        :return: A simplified dict_var based on the args
        """
        dict_var_simplified = {}

        for item in list_keys:
            dict_var_simplified[item] = dict_var[item]

        return dict_var_simplified

    def write_csv(self, path, tec, key):
        """
        Based on a dict (TEC), select the 'key' for printing all PRNs in columns in a CSV file format
        :param path: Path to save the CSV file
        :param tec: Dict object
        :param key: Key to be printed
        :return:
        """
        time = tec['time']

        f = csv.writer(open(path, "w+"))

        columns = []
        columns.append('time')
        for prn in tec[key]:
            columns.append(prn)

        f.writerow(columns)

        prns = list(tec[key].keys())
        elements = len(tec[key][prns[0]])
        for element in range(elements):
            row = []
            row.append(time[element])
            for prn in tec[key]:
                row.append(tec[key][prn][element])
            f.writerow(row)

    def process_tec_file(self, rinex_folder, file):
        """
        TEC workflow, with the calculus of the whole day. At the end of this workflow, the

        :param rinex_folder: Absolute folder with a set of rinex files
        :param file: the rinex name file
        :return: All the TEC estimative, including the piercing points, slant factor, daily TEC, bias receiver,
        detrended, relative, absolute and, finally, vertical TEC
        """
        tec = {}

        input_files = helper.InputFiles()
        utils = helper.Utils()
        cycle_slip = cs.CycleSlip()
        tec_estimation = calc.TECEstimation()
        bias_estimation = calc.BiasEstimation()
        quality_control = calc.QualityControl()

        start = time.perf_counter()

        logging.info("- {} - TEC by fractions of {} {} a day, and bias receiver estimation".
                     format(file, settings.TEC_RESOLUTION_VALUE, settings.TEC_RESOLUTION))
        logging.info("Preparing inputs...")
        hdr, obs, orbit, dcb, factor_r, l1_col, l2_or_l3_col, p1_or_c1_col, p2_or_c2_col, l2_channel, \
        constellations = input_files.prepare_inputs(rinex_folder, file)

        tec['metadata'] = {
            'creation-date': datetime.datetime.utcnow(),
            'modification-date': datetime.datetime.utcnow(),
            'rinex-file': file,
            'rinex-path': os.path.join(rinex_folder, file),
            'rinex-date': utils.rinex_str_date_to_datetime(hdr),
            'rinex-precision': hdr['interval'],
            'keys-saved': settings.KEYS_SAVE,
            'constellations-desired': settings.CONSTELATIONS,
            'constellations-used-in-the-calc': constellations,
            'min-elev-angle': settings.MIN_ELEVATION_ANGLE,
            'station': hdr['MARKER NAME'].strip(),
            'orbit-name': orbit['path'],
            'dcb-name': dcb['path']
        }

        logging.info("Converting timestamp str in datetimes...")
        tec['time'] = utils.array_timestamp_to_datetime(obs.time)

        try:
            logging.info("Calculating slant factor for DTEC calculation...")
            tec['slant-dtec'] = tec_estimation.slant(hdr, obs, orbit, 0)

            logging.info("Calculating slant factor for TEC estimation...")
            tec['slant'] = tec_estimation.slant(hdr, obs, orbit, 1)

            logging.info("Correcting Cycle-Slip...")
            tec['relative-l1-l2'] = cycle_slip.cycle_slip_analysis(obs, tec, factor_r, l1_col, l2_or_l3_col,
                                                                   p1_or_c1_col, p2_or_c2_col, l2_channel)

            logging.info("Calculating dTEC...")
            tec['detrended'] = tec_estimation.detrended(tec, factor_r, l2_channel)

            logging.info("Calculating relative TEC...")
            tec['relative'] = tec_estimation.relative(tec, obs, factor_r, dcb, p1_or_c1_col, p2_or_c2_col)

            logging.info("Estimating TEC and Bias with daily measurements...")
            tec['matrixes'], tec['bias'] = bias_estimation.estimate_bias(tec, constellations)

            logging.info("Quality control - calculating the quality of the estimates...")
            tec['quality'] = quality_control.quality_control(tec['matrixes'], tec['bias'])

            logging.info("Filtering - filtering measures and error detection...")
            # tec['bias'] = quality_control.check_quality(obs, tec, constellations, rinex_folder, file)

            logging.info("Calculating absolute TEC...")
            tec['absolute'] = tec_estimation.absolute(tec, constellations)

            logging.info("Calculating vertical TEC...")
            tec['vertical'] = tec_estimation.vertical(tec, orbit)

            # utils.write_csv(tec, 'absolute.csv', tec[absolute])

        except Exception as e:
            logging.error(':::: EXCEPTION thrown during {} processing: {}! File skipped!\n'.format(file, e))

        stop = time.process_time()
        logging.info("Processing done for {}! Time: {} minutes".format(file, float((start - stop) / 60)))
        logging.info("-----------------------------------------------------------------------------------")

        tec_to_be_stored = self.simplify_dict(tec, settings.KEYS_SAVE)

        return tec_to_be_stored

