"""
Commom settings to all applications
"""
import datetime

from decouple import config

PATH_DCB = config('PATH_DCB', default='/embrace/tec/dcb/')
PATH_ORBIT = config('PATH_ORBIT', default='/embrace/tec/orbit/')
PATH_GLONASS_CHANNEL = config('PATH_GLONASS_CHANNEL', default='/embrace/tec/glonasschannel/')

URL_DCB_MGEX = 'ftp://cddis.gsfc.nasa.gov/pub/gps/products/mgex/dcb/'
URL_ORBIT_MGEX = 'ftp://cddis.gsfc.nasa.gov/pub/gps/products/mgex/'
URL_GLONASS_CHANNELS = 'https://www.glonass-iac.ru/en/CUSGLONASS/'

VIEW_ANGLE = 30
INTERVAL_IN_SECS = 30
REQUIRED_VERSION = 3.01
TEC_RESOLUTION_ESTIMATION = datetime.timedelta(hours=1)

A = 40.3
TECU = 1.0e16
C = 299792458
DIFF_TEC_MAX = 0.05
LIMIT_STD = 7.5

FREQUENCIES = {'G': {'1': 1.57542e9, '2': 1.22760e9}, 'R': {'1': 1602.0e+6, '2': 1246.0e+6}}

THRESHOLD_VARIANCE_BIAS = 2
INITIAL_HOUR_RECALC_BIAS = 7
FINAL_HOUR_RECALC_BIAS = 17

GNSS_EPOCH = 3657

ELLIPTICITY = 298.257223563
RADIUS_EQUATOR = 6378.137

CONSTELATIONS = ['R']

ALT_IONO = 6670.0
ALT_IONO_BOTTOM = 6620.0
ALT_IONO_TOP = 6870.0

EARTH_RAY = 6371.0

COLUMNS_IN_RINEX = {'3.03': {'G': {'L1': 'L1C', 'L2': 'L2W', 'C1': 'C1C', 'P1': 'C1W', 'P2': 'C2W'},
                             'R': {'L1': 'L1C', 'L2': 'L2W', 'C1': 'C1C', 'P1': 'C1P', 'P2': 'C2P'}
                             },
                    '3.02': {'G': {'L1': 'L1', 'L2': 'L2', 'C1': 'C1C', 'P1': 'C1W', 'P2': 'C2W'},
                             'R': {'L1': 'L1', 'L2': 'L2', 'C1': 'C1C', 'P1': 'C1P', 'P2': 'C2P'}
                             },
                    '3.01': {'G': {'L1': 'L1', 'L2': 'L2', 'C1': 'C1C', 'P1': 'C1W', 'P2': 'C2W'},
                             'R': {'L1': 'L1', 'L2': 'L2', 'C1': 'C1C', 'P1': 'C1P', 'P2': 'C2P'}
                             }
                    }

OBS_MGEX = [
    ('C1-P1', (('G', ('C1C', 'C1W')), ('R', ('C1C', 'C1P')))),
    ('P1-P2', (('G', ('C1W', 'C2W')), ('R', ('C1P', 'C2P'))))
]

REGEX_DCB_1 = r'\s\w{3}\s[\s\w\d]{5}(\s\w\d\d)[\s\w\d]+'
REGEX_DCB_2 = r'\s(\s[\d]{2}:\d{3}:\d{5}){2}[\s\w]{4}([-\d\s]{18}.\d{4})([-\d\s]{7}.\d{4})'

REGEX_GLONASS_CHANNEL = r'\d\/([\s\d]{3})\|([-\s\d]{6})'
REGEX_GLONASS_CHANNEL_RINEX = r'R(\d\d)([-\d\s]{4})'
REGEX_RINEX_DATE = r'(\d{4})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+.0)'
