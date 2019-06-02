<img src="sw_logo.png" width="64">

# TEC and Bias receiver estimation
[![](https://img.shields.io/github/license/embrace-inpe/swds-api-downloader.svg)](https://github.com/embrace-inpe/swds-api-downloader/blob/master/LICENSE)
[![](https://img.shields.io/badge/python-3-blue.svg)](https://www.python.org/)
[![](https://img.shields.io/badge/INPE-EMBRACE-orange.svg)](http://www2.inpe.br/climaespacial/portal/pt/)
[![](https://img.shields.io/badge/coverave-20%25-orange.svg)](https://github.com/embrace-inpe/swds-api-downloader)

Application TEC (by station) and receiver bias estimation through GNSS data processing and analysis, used for monitoring of ionospheric terrestrial layer. 

Please, fell free to read more in [EMBRACE](http://www2.inpe.br/climaespacial/portal/pt/).

#### Contributors involved
###### Departamento de Ciências Espaciais e Atmosféricas (CEA-II) - INPE
Dr. Cristiano Max Wrasse (Pesquisador) [_cristiano.wrasse@inpe.br_]  
Dr. Cosme A. O. B. Figueiredo [_cosme.figueiredo@inpe.br_]  

###### Desenvolvimento - EMBRACE/INPE
Dr. Rodolfo G. Lotte [_rodolfo.lotte@inpe.br_]  
 
***
## GNSS's versions file covered:
- GNSS rinex version 3.01
- GNSS rinex version 3.02
- GNSS rinex version 3.03

PS:. The EMBRACE TEC and Bias estimation make use of external library `georinex` for reading rinex files, which includes 
other versions! To read more about, check the [link](https://pypi.org/project/georinex/).

## Features
This module includes the calculation of:

- Cycle-Slip correction
- Relative TEC
- Slant Factor
- Daily TEC and receiver bias estimation
- Absolute TEC
- Vertical TEC

## Output
This module runs for a set of rinex files. For each rinex processed, a 
python dictionary JSON-like structured is generated.

## Contributing:
### TODO list: 
Of course, the model is far completed. Besides the TODOs marks into the code, the analysis for another 
situations are still demanded. If you notice some kind of mistake in the code, or just notice that could improve or 
optimize any method, please, fell free to clone this project and help our team to get better estimates.

### Preparing your environment:
To contribute, you have to clone this repository and start your analysis/debugs/so on.

After to clone the repository [git](https://github.com/embrace-inpe/tec) in your local directory, 
the following the instructions will guide you to execute the model in your computer.

1. [Create a isolated Python environment with virtualenv](#2 )
2. [Installing dependencies with `pip`](#2-Instalando-dependncias)
3. [Setting up the constants](#4-Configurando-constantes)
4. [Executing the `main.py`](#5-Execuo-do-programa)
5. [Performance](#6-Performance)

#### 1. Create a isolated Python environment with virtualenv
Go to `tec/` and create an isolated environment with:

```console
$ python -m venv .venv
```

To activate your `venv`, type the following: 

```console
$ source .venv/bin/activate
```

If everything is ok, you will see the virtual env prefix:
```console
(.venv) usuario@maquina:<seu-diretorio>$ your commands here
```

Now, everything installed through `pip`, will be encapsulated by the `venv`, and will be available only in this 
respective scope.

To deactivate, just type: 

```console
$ deactivate
```

#### 2. Installing the dependencies with `pip`
After you have set all the variables in `.env` with you personal information, you need to install the 
dependencies listed in `requirements.txt`:
```
pip3 install -r requirements.txt
```
Then, run:
```
python3 runtests.py
```

#### 3. Setting up the constants in `.env`
The `settings.py` is the module responsible to store all the constants used by the model. Besides these constants, 
another ones are also essential to make the process possible, the ones respected to your 
computer!

So, in this case, the constants that change the computer by computer, is setup in `.env` file. 
Here, you can set what you want to consider in the modelling. For instance, the 
constellations, the minimum required rinex version, the resolution of estimation (1 hour, 
10 minutes, so on), and the paths the rinex are! 

First, copy the file `.env-example` and rename for `.env`. You will see all variables setted like:
```
PATH_DCB=
PATH_ORBIT=
PATH_GLONASS_CHANNEL=
RINEX_FOLDER=

MIN_REQUIRED_VERSION=2.11
CONSTELATIONS=['G', 'R', 'E', 'C']
TEC_RESOLUTION=hours
TEC_RESOLUTION_VALUE=1
KEYS_SAVE=['time', 'slant-dtec', 'slant', 'detrended', 'bias', 'quality', 'vertical']
```

- `PATH_DCB`:
- `PATH_ORBIT`: 
- `PATH_GLONASS_CHANNEL`: 
- `RINEX_FOLDER`: 
- `MIN_REQUIRED_VERSION`: 
- `TEC_RESOLUTION`: 
- `TEC_RESOLUTION_VALUE`: 
- `KEYS_SAVE`:
 
#### 4. Executing the `main.py`
Once the `.venv` activated, the main execution should be made by the `main.py` script:
```console
$ python main.py
```

If all correct, the output will be something like:
```console
[2019.06.02 13:34:01] {tec.py         :81  } INFO : - ango2220.14o - TEC by fractions of 1 hours a day, and bias receiver estimation 
[2019.06.02 13:34:01] {tec.py         :82  } INFO : Preparing inputs... 
[2019.06.02 13:34:01] {helper.py      :843 } INFO : >> Validating file and measures... 
[2019.06.02 13:34:01] {helper.py      :845 } INFO : >>>> Rinex version 3.01! 
[2019.06.02 13:34:01] {helper.py      :756 } INFO : >>>>>> Column L1 is not available for constellation E. TEC wont be consider for this constellation! 
[2019.06.02 13:34:01] {helper.py      :756 } INFO : >>>>>> Column L1 is not available for constellation C. TEC wont be consider for this constellation! 
[2019.06.02 13:34:01] {helper.py      :857 } INFO : >> Reading rinex measures... Only constellation(s) ['G'] will be considered! 
[2019.06.02 13:34:35] {helper.py      :882 } INFO : >> Downloading GLONASS channels in https://www.glonass-iac.ru/en/CUSGLONASS/ for pseudorange calculation... 
[2019.06.02 13:34:35] {downloads.py   :78  } INFO : >>>> Glonass Channel File already exist. Download skipped! 
[2019.06.02 13:34:35] {parser.py      :166 } INFO : >> Starting GLONASS channels parsing... 
[2019.06.02 13:34:35] {helper.py      :918 } INFO : >> Downloading Orbit files...
[...]
```

#### 5. Performance


## Log
Download errors will be listed in tec.log on root path.

## Help
Any question or suggestion, please, contact our support sending an email to `desenvolvimento.emabrace@inpe.br` or any 
of the contributers.


