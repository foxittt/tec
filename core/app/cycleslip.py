import datetime
import logging

import numpy as np
import settings as settings

from app.utils.helpers import Utils
from scipy.signal import find_peaks


class CycleSlip:
    """
    If there is some cycle-slip loses, it will return a Python obs object with the cycles preserved
    """
    requiried_version = str(settings.REQUIRED_VERSION)
    cols_var = settings.COLUMNS_IN_RINEX[requiried_version]

    def _detect(self, rtec_no_nan, prn):
        """
        The method precisely detect the variations over relative TEC (by the L1 and L2 differences). When a "degree"
        effect is present, a peak on the fourth derivative occurs. The index of this peak is store on "indexes"
        variable and then returned

        :param rtec: The relative TEC (with NaNs)
        :return: The rtec indexes (in the array with NaN values) where the "degree" effect is presented
        """
        fourth_der_not_nan = np.diff(rtec_no_nan, n=4)
        std_fourth_der_not_nan = np.nanstd(fourth_der_not_nan) * settings.LIMIT_STD
        indexes = find_peaks(abs(fourth_der_not_nan), height=std_fourth_der_not_nan)[0]
        indexes = np.array(indexes)

        if len(indexes) == 0:
            logging.info(">>>>>> No discontinuities detected (by final differences) for PRN {}".format(prn))
        else:
            logging.info(">>>>>> Discontinuities detected in {} (not NaN) for PRN {}".format(indexes, prn))

        return indexes

    def _correct(self, l1, l2, c1, p2, rtec, mwlc, f1, f2, factor_1, factor_2, index):
        """
        As a result of the detection of inconsistencies, this method compensate the irregularities on L1 and L2
        observations, directly at the rinex (obs variable)

        :param l1: L1 measures (no NaN values)
        :param l2: L2 measures (no NaN values)
        :param c1: C1 measures (no NaN values)
        :param p2: P2 measures (no NaN values)
        :param rtec_not_nan: relative TEC (no NaNs)
        :param mlwc: No NaN MLWC factor
        :param f1: F1 frequency (either GPS or GLONASS)
        :param f2: F2 frequency (either GPS or GLONASS)
        :param factor_1: first factor of calculus (either GPS or GLONASS)
        :param factor_2: second factor of calculus (either GPS or GLONASS)
        :param index: The point of inconsistence (array with no NaN values)
        :return: The corrected observation file and the respective relative TEC
        """
        diff_rtec = rtec[index] - rtec[index - 1]
        diff_mwlc = mwlc[index] - mwlc[index - 1]

        var_1 = diff_mwlc * settings.C
        var_2 = var_1 / f1
        diff_2 = round((diff_rtec - var_2) * factor_2)
        diff_1 = diff_2 + round(diff_mwlc)

        cor_r_1 = l1[index] - diff_1
        cor_r_2 = l2[index] - diff_2

        var_corr_1 = cor_r_1 / f1
        var_corr_2 = cor_r_2 / f2

        rtec[index] = (var_corr_1 - var_corr_2) * settings.C
        mwlc[index] = (cor_r_1 - cor_r_2) - (f1 * c1[index] + f2 * p2[index]) * factor_1

        l1[index:len(l1)] -= diff_1
        l2[index:len(l2)] -= diff_2

        return rtec, l1, l2

    def _detect_and_correct_cycle_slip(self, obs_time, l1, l2, c1, p2, f1, f2, factor_1, factor_2, prn):
        """
        Start the variables to check cycle-slip for each PRN presented in rinex files

        :param obs_time: The array of all times (already converted to datetimes) regarding the current rinex file
        :param l1: L1 measures (with NaN values)
        :param l2: L2 measures (with NaN values)
        :param c1: C1 measures (with NaN values)
        :param p2: P2 measures (with NaN values)
        :param f1: F1 frequency (either GPS or GLONASS)
        :param f2: F2 frequency (either GPS or GLONASS)
        :param factor_1: first factor of calculus (either GPS or GLONASS)
        :param factor_2: second factor of calculus (either GPS or GLONASS)
        :param prn: The respective PRN
        :return: The relative TEC base on the differences between L1 and L2 (rtec_nan), with cycle-slip corrections
        """
        j_start = 0

        rtec_nan = ((l1 / f1) - (l2 / f2)) * settings.C
        mwlc_nan = (l1 - l2) - (f1 * c1 + f2 * p2) * factor_1

        rtec_no_nan = rtec_nan.copy()
        mwlc_no_nan = mwlc_nan.copy()
        rtec_no_nan = rtec_no_nan[~np.isnan(rtec_no_nan)]
        mwlc_no_nan = mwlc_no_nan[~np.isnan(mwlc_no_nan)]

        nan_pos = np.where(np.isnan(rtec_nan))
        nan_pos = np.array(nan_pos).flatten().tolist()
        not_nan_pos = np.where(~np.isnan(rtec_nan))
        not_nan_pos = np.array(not_nan_pos).flatten().tolist()

        not_nan_pos = np.array(not_nan_pos).flatten().tolist()
        not_nan_obs_time = [obs_time[x] for x in not_nan_pos]
        l1_not_nan = [l1[x] for x in not_nan_pos]
        l2_not_nan = [l2[x] for x in not_nan_pos]
        c1_not_nan = [c1[x] for x in not_nan_pos]
        p2_not_nan = [p2[x] for x in not_nan_pos]

        logging.info(">>>> Detecting peaks on the 4th order final differences in rTEC...")
        indexes = self._detect(rtec_no_nan, prn)

        logging.info(">>>> Finding discontinuities and correcting cycle-slips (PRN {})...".format(prn))
        for i in range(1, len(not_nan_obs_time)):
            rtec_no_nan[i] = ((l1_not_nan[i] / f1) - (l2_not_nan[i] / f2)) * settings.C
            mwlc_no_nan[i] = (l1_not_nan[i] - l2_not_nan[i]) - (f1 * c1_not_nan[i] + f2 * p2_not_nan[i]) * factor_1

            t1 = not_nan_obs_time[i - 1]
            t2 = not_nan_obs_time[i]

            if t2 - t1 > datetime.timedelta(minutes=15):
                j_start = i
                continue

            if i in indexes:
                logging.info(">>>>>> Indexes match ({}): correcting cycle-slips...".format(i, prn))
                rtec_no_nan, l1_not_nan, l2_not_nan = self._correct(l1_not_nan, l2_not_nan, c1_not_nan, p2_not_nan,
                                                                    rtec_no_nan, mwlc_no_nan, f1, f2,
                                                                    factor_1, factor_2, i)

            if i - j_start + 1 >= 12:
                add_tec = 0
                add_tec_2 = 0

                for jj in range(1, 10):
                    add_tec = add_tec + rtec_no_nan[i - jj] - rtec_no_nan[i - jj - 1]
                    add_tec_2 = add_tec_2 + pow(rtec_no_nan[i - jj] - rtec_no_nan[i - jj - 1], 2)

                p_mean = add_tec / 10
                p_dev = np.maximum(np.sqrt(add_tec_2 / 10 - pow(p_mean, 2)), settings.DIFF_TEC_MAX)
            else:
                p_mean = 0
                p_dev = settings.DIFF_TEC_MAX * 2.5

            pmin_tec = p_mean - p_dev * 2
            pmax_tec = p_mean + p_dev * 2
            diff_rtec = rtec_no_nan[i] - rtec_no_nan[i - 1]

            if not pmin_tec < diff_rtec and diff_rtec <= pmax_tec:
                rtec_no_nan, l1_not_nan, l2_not_nan = self._correct(l1_not_nan, l2_not_nan, c1_not_nan, p2_not_nan,
                                                                    rtec_no_nan, mwlc_no_nan, f1, f2,
                                                                    factor_1, factor_2, i)

        # np_zeros = np.zeros(len(nan_pos))
        # rtec_no_nan = np.concatenate((rtec_no_nan, np_zeros), axis=0)
        # rtec_no_nan[nan_pos] = np.nan

        return rtec_no_nan

    def cycle_slip_analysis(self, obs, tec, factor_glonass):
        """
        Correct the lost of signals during the communications between satellites and receiver, usually called cycle-slip
        correction. It improves and help to get a more consistent estimate

        :param obs:
        :param tec:
        :param factor_glonass:
        :return:
        """
        obs_time = tec['time']
        slant = tec['slant']
        prns = obs.sv.values

        for prn in prns:
            # TODO: reduzir janela de correção
            # slant_reduced = slant[prn][2] > 30

            l1 = np.array(obs[self.cols_var[prn[0:1]]['L1']].sel(sv=prn).values)
            l2 = np.array(obs[self.cols_var[prn[0:1]]['L2']].sel(sv=prn).values)
            c1 = np.array(obs[self.cols_var[prn[0:1]]['C1']].sel(sv=prn).values)
            p2 = np.array(obs[self.cols_var[prn[0:1]]['P2']].sel(sv=prn).values)

            if prn[0:1] == 'R':
                f1 = factor_glonass[prn[1:]][0]
                f2 = factor_glonass[prn[1:]][1]
                factor_1 = factor_glonass[prn[1:]][2]
                factor_2 = factor_glonass[prn[1:]][3]
            elif prn[0:1] == 'G':
                f1 = settings.F1
                f2 = settings.F2
                factor_1 = settings.factor_1
                factor_2 = settings.factor_2

            rtec_corrected = self._detect_and_correct_cycle_slip(obs_time, l1, l2, c1, p2,
                                                                 f1, f2, factor_1, factor_2, prn)

        return obs
