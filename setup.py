
from setuptools import setup

setup(version='0.0.2',
      name='mendelianerror',
      py_modules=['mendelianerror'],
      description="probability of mendelian error in trios",
      entry_points={
      'console_scripts': ['mendelianerror = mendelianerror:_main']},
      long_description=open('README.md').read(),
      author="Brent Pedersen",
      author_email="bpederse@gmail.com",
      zip_safe=False,
      classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: MIT License',
              'Topic :: Scientific/Engineering :: Bio-Informatics'])
