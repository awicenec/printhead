from setuptools import setup, find_packages

setup(name='printhead',
      version='5.0',
      packages=find_packages(),
      install_requires=[],
      entry_points={'console_scripts': [
          'printhead = printhead.__main__:main',]
      })