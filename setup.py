from setuptools import setup, find_packages

setup(name='tec-embrace',
      version='0.1',
      url='https://github.com/embrace-inpe/tec',
      license='MIT',
      author='Estudo e Monitoramento Brasileiro do Clima Espacial (EMBRACE/INPE)',
      author_email='desenvolvimento.embrace@gmail.com',
      description='Add static script_dir() method to Path',
      packages=find_packages(exclude=['tests']),
      long_description=open('README.md').read(),
      zip_safe=False)