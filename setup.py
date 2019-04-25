from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name = 'mecan4cna',
     version = '0.17',
     description = 'Minimum Error Calibration and Normalization for Copy Number Analysis',
     long_description = readme(),
     license = 'MIT',
     url='https://github.com/baudisgroup/segment-liftover',
     author = 'Bo Gao',
     author_email = 'kaye_gao@hotmail.com',
     classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    keywords = 'calibration normalization copy number',
    packages = ['mecan4cna'],
    install_requires = [
        'click',
        'pandas',
        'numpy',
        'matplotlib'
        ],
    python_requires = '>=3.6',
    entry_points = {
        'console_scripts': ['mecan4cna = mecan4cna.cli:main']
    },
    include_package_data = True
)