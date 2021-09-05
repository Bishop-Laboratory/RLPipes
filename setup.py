from setuptools import setup, find_packages

setup(
    name='rlpipes',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'pysradb',
        'pysam',
        'pyfastx',
        'pandas==1.2.0'
    ],
    entry_points={
        'console_scripts': [
            'RLPipes = rlpipes.cli:cli',
        ],
    }
)
