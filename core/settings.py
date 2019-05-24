"""
Commom settings to all applications
"""
import datetime

from decouple import config

RINEX_FOLDER = config('RINEX_FOLDER', default='/embrace/cosme/rinex/301/ango_2014/')
PATH_DCB = config('PATH_DCB', default='/embrace/tec/dcb/')
PATH_ORBIT = config('PATH_ORBIT', default='/embrace/tec/orbit/')
PATH_GLONASS_CHANNEL = config('PATH_GLONASS_CHANNEL', default='/embrace/tec/glonasschannel/')

URL_DCB_MGEX = 'ftp://cddis.gsfc.nasa.gov/pub/gps/products/mgex/dcb/'
URL_ORBIT_MGEX = 'ftp://cddis.gsfc.nasa.gov/pub/gps/products/mgex/'
URL_GLONASS_CHANNELS = 'https://www.glonass-iac.ru/en/CUSGLONASS/'

VIEW_ANGLE = 30
INTERVAL_IN_SECS = 30
MIN_REQUIRED_VERSION = 3.01
TEC_RESOLUTION_ESTIMATION = datetime.timedelta(hours=1)

KEYS_SAVE = ['time', 'slant', 'detrended', 'bias', 'quality', 'vertical']
CONSTELATIONS = ['G', 'E', 'C']
MIN_ELEVATION_ANGLE = 30
A = 40.3
TECU = 1.0e16
C = 299792458
DIFF_TEC_MAX = 0.05
LIMIT_STD = 7.5

FREQUENCIES = {'G': [1.57542e9, 1.22760e9, 1.17645e9],
               'R': [1602.0e+6, 1246.0e+6, 1.20202e9],
               'E': [1.57542e9, 1.17645e9, 1.20714e9],
               'C': [1.56109e9, 1.20714e9, 1.26852e9]}

THRESHOLD_VARIANCE_BIAS = 2
INITIAL_HOUR_RECALC_BIAS = 7
FINAL_HOUR_RECALC_BIAS = 17

GNSS_EPOCH = 3657
ELLIPTICITY = 298.257223563
RADIUS_EQUATOR = 6378.137
ALT_IONO = 6670.0
ALT_IONO_BOTTOM = 6620.0
ALT_IONO_TOP = 6870.0
EARTH_RAY = 6371.0

COLUMNS_IN_RINEX = {'3.03': {'G': {'L1': 'L1C', 'L2': 'L2W', 'L3': 'L5Q',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1W', 'P2': 'C2W'},
                             'R': {'L1': 'L1C', 'L2': 'L2C', 'L3': 'L3I',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'},
                             'E': {'L1': 'L1C', 'L2': 'L5Q', 'L3': 'L7Q',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'},
                             'C': {'L1': 'L2I', 'L2': 'L7I', 'L3': 'L6I',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'}
                             },
                    '3.02': {'G': {'L1': 'L1', 'L2': 'L2', 'L3': 'L2N',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1W', 'P2': 'C2W'},
                             'R': {'L1': 'L1', 'L2': 'L2', 'L3': 'L3I',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'},
                             'E': {'L1': 'L1C', 'L2': 'L2W', 'L3': 'L7I',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'},
                             'C': {'L1': 'L1C', 'L2': 'L2W', 'L3': 'L6Q',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'}
                             },
                    '3.01': {'G': {'L1': 'L1', 'L2': 'L2', 'L3': 'L5Q',
                                   'C1': 'C1C', 'C2': 'C2W', 'P1': 'C1W', 'P2': 'C2W'},
                             'R': {'L1': 'L1', 'L2': 'L2', 'L3': 'XXX',
                                   'C1': 'C1C', 'C2': 'C2P', 'P1': 'C1P', 'P2': 'C2P'},
                             'E': {'L1': 'L1', 'L2': 'L2', 'L3': 'L7I',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'},
                             'C': {'L1': 'L1', 'L2': 'L2', 'L3': 'XXX',
                                   'C1': 'C1C', 'C2': 'C2C', 'P1': 'C1P', 'P2': 'C2P'}
                             }
                    }

# TODO: Completar colunas de E. MGEX também obedece às diferentes versões? Se sim, incluir como em COLUMNS_IN_RINEX
OBS_MGEX = [
    ('C1-P1', (('G', ('C1C', 'C1W')), ('R', ('C1C', 'C1P')), ('E', ('C1C', 'C5Q')), ('C', ('C2I', 'C7I')))),
    ('P1-P2', (('G', ('C1W', 'C2W')), ('R', ('C1P', 'C2P')))),
    ('C2-P2', (('G', ('C2C', 'C2W')), ('R', ('C2C', 'C2P'))))
]

REGEX_DCB_1 = r'\s\w{3}\s[\s\w\d]{5}(\s\w\d\d)\s+'
REGEX_DCB_2 = r'\s+(\s[\d]{2,4}:\d{3}:\d{5}){2}[\s\w]{4}([-\d\s]{18}.\d{4})([-\d\s]{7}.\d{4})'

REGEX_GLONASS_CHANNEL = r'\d\/([\s\d]{3})\|([-\s\d]{6})'
REGEX_GLONASS_CHANNEL_RINEX = r'R(\d\d)([-\d\s]{4})'
REGEX_RINEX_DATE = r'(\d{4})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+.0)'
