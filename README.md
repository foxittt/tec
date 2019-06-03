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

###### Development team - EMBRACE/INPE
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
This module runs for a set of rinex files. For each rinex processed, a python dictionary JSON-like structured is generated.

## Contributing:
### TODO list: 
Of course, the model is far completed. Besides the TODOs marks into the code, the analysis for another 
situations are still demanded. If you notice some kind of mistake in the code, or just notice that could improve or 
optimize any method, please, fell free to clone this project and help our team to get better estimates.

### Preparing your environment:
To contribute, you have to clone this repository and start your analysis/debugs/ so on.

After to clone the repository [git](https://github.com/embrace-inpe/tec) in your local directory, 
the following the instructions will guide you to execute the model in your computer.

#### 1. Updating your Python
The EMBRACE TEC and Bias estimation was developed using version Python 3.7+. If you do not have this version in your computer, an upgrade will be need. 

First, check the version of your python typing `python -V`. If old versions are diplayed, thus, follow the steps below. For most of Ubuntu users, two versions are available on the OS, `python2.7`, and `python3`, which still do not fill the requirement due the `python3` be usually under version 3.6 (to check this, please, type `python3 -V`). If not, the following steps can solve the problem.

To start the upgrade, update your system:
```
sudo apt-get update && sudo apt-get upgrade && sudo apt-get update
```
then, run:
```
sudo apt-get install python3.7
```

If everything is ok, edit your `bashrc` file:
```
sudo nano ~/.bashrc
```
and update the file:
```
source ~/.bashrc
```

and include at the last line: `alias python=python3.7`. Type:
```
python -V
``` 

## Example of use
An example of using this package/library, can be found (here)[https://github.com/rodolfolotte/tec-example].

## Log
Download errors will be listed in tec.log on root path.

## Help
Any question or suggestion, please, contact our support sending an email to `desenvolvimento.emabrace@inpe.br` or any 
of the contributers.


