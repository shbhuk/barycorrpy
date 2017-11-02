from setuptools import setup

setup(name='barycorrpy',
      version='0.1',
      description='Barycentric Velocity correction at 1 cm/s level',
      url='https://github.com/shbhuk/barycorrpy',
      author='Shubham Kanodia',
      author_email='shbhuk@gmail.com',
      install_requires=['astropy','jplephem'],
      classifiers=['Topic :: Scientific/Engineering :: Astronomy'],
      include_package_data=True
      )
      
