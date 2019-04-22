import argparse
import logging
import os
import sys
import time
from logging import handlers

import settings as settings
import helper as helper
import estimate as calc
import app.cycleslip as cs


def tec_by_station(rinex_folder):
    """
    Main workflow of TECMAP by station. This routine can sequentially calculate the TEC relative, the slant factor,
    absolute, and vertical. Besides, it can generate an accurate estimate of bias receptor along the day

    :param rinex_folder: Absolute path to the rinex's folder (files respecting to N-stations)
    :return: None
    """
    tec = {}
    input_files = helper.InputFiles()
    utils = helper.Utils()
    cycle_slip = cs.CycleSlip()
    tec_estimation = calc.TECEstimation()
    bias_estimation = calc.BiasEstimation()
    quality_control = calc.QualityControl()

    files = sorted(os.listdir(rinex_folder))

    for file in files:
        start = time.perf_counter()

        logging.info("- {} - TEC by fractions of {} a day, and bias receiver "
                     "estimation -----------".format(file, settings.TEC_RESOLUTION_ESTIMATION))
        logging.info("Preparing inputs:")
        hdr, obs, orbit, dcb, factor_r = input_files.prepare_inputs(rinex_folder, file)

        tec['time'] = utils.array_timestamp_to_datetime(obs.time)

        logging.info("Calculating slant factor:")
        tec['slant'] = tec_estimation.slant(hdr, obs, orbit)

        logging.info("Correcting Cycle-Slip...")
        tec['relative-l1l2'] = cycle_slip.cycle_slip_analysis(obs, tec, factor_r)

        logging.info("Calculating relative TEC:")
        tec['relative-p1c1'] = tec_estimation.relative(hdr, obs, factor_r, dcb)

        logging.info("Estimating TEC and Bias with daily measurements:")
        bias = bias_estimation.estimate_bias(tec)

        logging.info("Quality control - calculating the quality of the estimates:")
        tec['quality'] = quality_control.quality_control(bias)

        logging.info("Filtering - filtering measures and error detection:")
        tec['reestimated'] = quality_control.check_quality(tec, rinex_folder, file)

        logging.info("Calculating absolute TEC:")
        tec['absolute'] = tec_estimation.absolute(tec, bias)

        logging.info("Calculating vertical TEC:")
        tec['vertical'] = tec_estimation.vertical(tec)

        stop = time.process_time()
        logging.info("Processing done for {}! Time: {} minutes".format(file, float((start - stop) / 60)))
        logging.info("-----------------------------------------------------------------------------------------------")


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='EMBRACE/INPE TEC map and receptor bias estimation')

    parser.add_argument('-rinex_folder', action="store", dest='rinex_folder',
                        help='Rinex folder: formats 2.11, 3.03 and Hatanaka are accept.')
    parser.add_argument('-verbose', action="store", dest='verbose', help='Print log of processing.')

    args = parser.parse_args()

    if args.verbose:
        log = logging.getLogger('')
        log.setLevel(logging.INFO)
        format = logging.Formatter("[%(asctime)s] {%(filename)-15s:%(lineno)-4s} %(levelname)-5s: %(message)s ",
                                   datefmt='%Y.%m.%d %H:%M:%S')

        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(format)

        fh = logging.handlers.RotatingFileHandler(filename='tecmap.log', maxBytes=(1048576 * 5), backupCount=7)
        fh.setFormatter(format)

        log.addHandler(ch)
        log.addHandler(fh)
    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    tec_by_station(args.rinex_folder)

