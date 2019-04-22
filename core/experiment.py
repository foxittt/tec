import pandas as pd 
import numpy as np
import logging
import re
import estimate
import datetime
import bisect

A = 40.3
TECU = 1.0e16

TEC_RESOLUTION_ESTIMATION = datetime.timedelta(minutes=10)

REGEX_GLONASS_CHANNEL = r'\d\/([\s\d]{3})\|([-\s\d]{6})'
REGEX_GLONASS_CHANNEL_RINEX = r'R(\d\d)([-\d\s]{4})'

OBS_MGEX = {'C1-P1': {'G': ['C1C', 'C1W'], 'R': ['C1C', 'C1P']},
            'P1-P2': {'G': ['C1W', 'C2W'], 'R': ['C1P', 'C2P']}}

REGEX_DCB_1 = r'\s\w{3}\s[\s\w\d]{5}(\s\w\d\d)[\s\w\d]{11}'
REGEX_DCB_2 = r'\s([\d\s]{5}:\d{3}:\d{5}){2}[\s\w]{4}([-\d\s]{18}.\d{4})([-\d\s]{7}.\d{4})'


def vertical(tec):
    """
    When calculated, the TEC is function of satellites incident angles, sometimes, at the horizon. The vertical
    TEC it is the process to remove this influence, bringing the TEC perpendicular to the receiver, called vertical
    TEC -> Vertical = Absolute / Slant

    :param tec: Dict with TEC python object
    :return: The updated TEC object, now, with vertical TEC calculated
    """
    tec_v = {}
    tec_s_short = {}

    time_in_rinex = tec['time']
    time_in_orbit = tec['slant']['G01'][0]

    max_date_in_rinex = max(d for d in time_in_rinex if isinstance(d, datetime.date))
    min_date_in_rinex = min(d for d in time_in_rinex if isinstance(d, datetime.date))

    lower = bisect.bisect_left(time_in_orbit, min_date_in_rinex)
    upper = bisect.bisect_right(time_in_orbit, max_date_in_rinex)

    for prn in tec['slant']:
        tec_s_short[prn] = [time_in_rinex, tec['slant'][prn][1][lower:upper]]

    for prn, values in tec['absolute'].items():
        absolute_np = np.array(tec['absolute'][prn])
        slant_np = np.array(tec_s_short[prn][1])
        tec_v[prn] = absolute_np / slant_np

    return tec_v


def _indexes(array_time):
    indexes_fraction = []
    indexes_fraction_aux = []

    fraction_limit = array_time[0] + TEC_RESOLUTION_ESTIMATION

    for i, item in enumerate(array_time):
        if item < fraction_limit:
            indexes_fraction_aux.append(i)
        else:
            fraction_limit = item + TEC_RESOLUTION_ESTIMATION
            indexes_fraction.append(indexes_fraction_aux)
            indexes_fraction_aux = []
            indexes_fraction_aux.append(i)

    indexes_fraction.append(indexes_fraction_aux)

    return indexes_fraction


def _build_coefficients(tec):
    """
    Build the coefficients of the equation system. These terms are defined through a Least-Square Fitting method,
    which consist in the minimization of set of unknown variable, dispose in a set of equations, so called,
    equation system (see Otsuka et al. A new Technique for mapping of TEC using GPS network in Japan). The coefficients,
    are values organized by hours, each hour will receive a mean value of a specific PRN. For instance, for hour '0h' of
    group_1, will receive an array with TOTAL_OF_SATELLITES positions, each position corresponds to a mean value of
    1 / tec slant, for group_2, the each position corresponds to a mean value of relative / tec slant

    :param tec: The TEC object, with relative and slant factor, calculated by PRN
    :return: The group 1 and 2 of coefficients, which is 1 / slant_factor, and
    tec_relative / slant_factor, respectively
        For example:
            coefficients = {
                        'group_1':
                                    {
                                    '0h': [G01_mean, G02_mean, ..., N_sat_mean],
                                    '1h': [G01_mean, G02_mean, ..., N_sat_mean],
                                    ...
                                    'INTERVAL_A_DAY': [G01_mean, G02_mean, ..., N_sat_mean],
                                    },
                        'group_2':
                                    {
                                    '0h': [G01_mean, G02_mean, ..., N_sat_mean],
                                    '1h': [G01_mean, G02_mean, ..., N_sat_mean],
                                    ...
                                    'INTERVAL_A_DAY': [G01_mean, G02_mean, ..., N_sat_mean],
                                    }
                            }
    """
    coefficients = {}

    group_1 = {}
    group_2 = {}

    indexes = _indexes(tec['time'])

    for i, ind in enumerate(indexes):
        group_1_aux = []
        group_2_aux = []

        print(i, ind)

        for prn in tec['slant']:
            elements_slant = np.take(tec['slant'][prn], ind)
            elements_relat = np.take(tec['relative'][prn], ind)
            _1_slant = np.divide(1, elements_slant)
            _relative_slant = np.divide(elements_relat, elements_slant)

            avg_sla = np.mean(_1_slant, dtype=np.float32)
            avg_rel = np.mean(_relative_slant, dtype=np.float32)
            group_1_aux.append(avg_sla)
            group_2_aux.append(avg_rel)

        group_1["every_" + str(TEC_RESOLUTION_ESTIMATION) + "_frac_" + str(i)] = group_1_aux
        group_2["every_" + str(TEC_RESOLUTION_ESTIMATION) + "_frac_" + str(i)] = group_2_aux

    coefficients['group_1'] = group_1
    coefficients['group_2'] = group_2

    return coefficients


if __name__ == '__main__':
    # myArray = ['G01', 'G02', 'G03', 'G05', 'R01', 'R02', 'R03', 'R04', 'R05']
    # myData = []
    # for x in range(3):
    #     myData.append(myArray)
    #
    # # permite o numpy imprimir no terminal as matrizes completas, sem abreviar
    # desired_width = 256
    # np.set_printoptions(linewidth=desired_width)
    # np.set_printoptions(threshold=np.inf)
    #
    # estimate_tools = estimate.BiasEstimation()
    #
    # diagonal_matrix = estimate_tools.build_diagonal_matrix(myArray)
    # f_matrix = estimate_tools._build_matrix_f(myData)
    # p_matrix = estimate_tools._build_matrix_p(myData)
    #
    # print(f'\nSample input data:\n {myArray}\n\n'
    #       f'List of lists:\n {myData}\n\n'
    #       f'Simple diagonal matrix from input data:\n\n{diagonal_matrix}\n\n'
    #       f'F Matrix built from a list of lists:\n\n{f_matrix}\n\n'
    #       f'P Matrix built from a list of lists:\n\n{p_matrix}\n\n')

    coefficients = {
        'group_1':
        {
            '0h': [0.7, 0.7, 0.7, 0.7, 0.7],
            '1h': [1.2, 1.2, 1.2, 1.2, 1.2],
            '3h': [3.2, 3.2, 3.2, 3.2, 3.2],
        },
        'group_2':
        {
            '0h': [0.7, 0.7, 0.7, 0.7, 0.7],
            '1h': [1.2, 1.2, 1.2, 1.2, 1.2],
            '3h': [3.2, 3.2, 3.2, 3.2, 3.2],
        }
    }

    tec = {}

    tec_time = [datetime.datetime(2016, 10, 5, 0, 0, 0), datetime.datetime(2016, 10, 5, 0, 1, 0),
                datetime.datetime(2016, 10, 5, 0, 2, 0), datetime.datetime(2016, 10, 5, 1, 0, 0),
                datetime.datetime(2016, 10, 5, 1, 1, 0), datetime.datetime(2016, 10, 5, 1, 2, 0),
                datetime.datetime(2016, 10, 5, 2, 0, 0), datetime.datetime(2016, 10, 5, 2, 1, 0),
                datetime.datetime(2016, 10, 5, 2, 2, 0)]

    tec_r = {'G01': [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
             'G02': [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
             'G03': [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3],
             'G04': [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]}

    tec_a = {'G01': [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
             'G02': [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
             'G03': [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3],
             'G04': [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]}

    tec_s = {'G01': [1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1],
             'G02': [1.2, 1.2, 1.2, 6.5, 5.2, 9.7, 1.2, 1.2, 1.2],
             'G03': [1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3],
             'G04': [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4]}

    # tec_s = {'G01':
    #             [
    #              [datetime.datetime(2016, 10, 5, 0, 0, 0), datetime.datetime(2016, 10, 5, 0, 1, 0),
    #              datetime.datetime(2016, 10, 5, 0, 2, 0), datetime.datetime(2016, 10, 5, 1, 0, 0),
    #              datetime.datetime(2016, 10, 5, 1, 1, 0), datetime.datetime(2016, 10, 5, 1, 2, 0),
    #              datetime.datetime(2016, 10, 5, 2, 0, 0), datetime.datetime(2016, 10, 5, 2, 1, 0),
    #              datetime.datetime(2016, 10, 5, 2, 2, 0), datetime.datetime(2016, 10, 5, 3, 0, 0),
    #              datetime.datetime(2016, 10, 5, 3, 0, 0), datetime.datetime(2016, 10, 5, 3, 1, 0),
    #              datetime.datetime(2016, 10, 5, 3, 2, 0), datetime.datetime(2016, 10, 5, 4, 0, 0),
    #              datetime.datetime(2016, 10, 5, 4, 1, 0), datetime.datetime(2016, 10, 5, 4, 2, 0)],
    #              [float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1]
    #             ],
    #          'G02': [
    #              [datetime.datetime(2016, 10, 5, 0, 0, 0), datetime.datetime(2016, 10, 5, 0, 1, 0),
    #               datetime.datetime(2016, 10, 5, 0, 2, 0), datetime.datetime(2016, 10, 5, 1, 0, 0),
    #               datetime.datetime(2016, 10, 5, 1, 1, 0), datetime.datetime(2016, 10, 5, 1, 2, 0),
    #               datetime.datetime(2016, 10, 5, 2, 0, 0), datetime.datetime(2016, 10, 5, 2, 1, 0),
    #               datetime.datetime(2016, 10, 5, 2, 2, 0), datetime.datetime(2016, 10, 5, 3, 0, 0),
    #               datetime.datetime(2016, 10, 5, 3, 0, 0), datetime.datetime(2016, 10, 5, 3, 1, 0),
    #               datetime.datetime(2016, 10, 5, 3, 2, 0), datetime.datetime(2016, 10, 5, 4, 0, 0),
    #               datetime.datetime(2016, 10, 5, 4, 1, 0), datetime.datetime(2016, 10, 5, 4, 2, 0)],
    #              [1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')]
    #             ],
    #          'G03':
    #              [
    #              [datetime.datetime(2016, 10, 5, 0, 0, 0), datetime.datetime(2016, 10, 5, 0, 1, 0),
    #               datetime.datetime(2016, 10, 5, 0, 2, 0), datetime.datetime(2016, 10, 5, 1, 0, 0),
    #               datetime.datetime(2016, 10, 5, 1, 1, 0), datetime.datetime(2016, 10, 5, 1, 2, 0),
    #               datetime.datetime(2016, 10, 5, 2, 0, 0), datetime.datetime(2016, 10, 5, 2, 1, 0),
    #               datetime.datetime(2016, 10, 5, 2, 2, 0), datetime.datetime(2016, 10, 5, 3, 0, 0),
    #               datetime.datetime(2016, 10, 5, 3, 0, 0), datetime.datetime(2016, 10, 5, 3, 1, 0),
    #               datetime.datetime(2016, 10, 5, 3, 2, 0), datetime.datetime(2016, 10, 5, 4, 0, 0),
    #               datetime.datetime(2016, 10, 5, 4, 1, 0), datetime.datetime(2016, 10, 5, 4, 2, 0)],
    #              [float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3]
    #             ],
    #          'G04':
    #              [
    #                  [datetime.datetime(2016, 10, 5, 0, 0, 0), datetime.datetime(2016, 10, 5, 0, 1, 0),
    #                   datetime.datetime(2016, 10, 5, 0, 2, 0), datetime.datetime(2016, 10, 5, 1, 0, 0),
    #                   datetime.datetime(2016, 10, 5, 1, 1, 0), datetime.datetime(2016, 10, 5, 1, 2, 0),
    #                   datetime.datetime(2016, 10, 5, 2, 0, 0), datetime.datetime(2016, 10, 5, 2, 1, 0),
    #                   datetime.datetime(2016, 10, 5, 2, 2, 0), datetime.datetime(2016, 10, 5, 3, 0, 0),
    #                   datetime.datetime(2016, 10, 5, 3, 0, 0), datetime.datetime(2016, 10, 5, 3, 1, 0),
    #                   datetime.datetime(2016, 10, 5, 3, 2, 0), datetime.datetime(2016, 10, 5, 4, 0, 0),
    #                   datetime.datetime(2016, 10, 5, 4, 1, 0), datetime.datetime(2016, 10, 5, 4, 2, 0)],
    #              [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')]
    #             ]
    #     }

    tec['time'] = tec_time
    tec['relative'] = tec_r
    tec['absolute'] = tec_a
    tec['slant'] = tec_s

    coeff = _build_coefficients(tec)



