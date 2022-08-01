import pathlib
from setuptools import setup, find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='rlpipes',
    version='0.9.2',
    description="A standardized R-loop-mapping pipeline",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/Bishop-Laboratory/RLPipes",
    author="Henry Miller",
    author_email="millerh1@livemail.uthscsa.edu",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9"
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'pysradb',
        'pysam>=0.17.0',
        'pyfastx',
        'pandas==1.2.0'
    ],
    entry_points={
        'console_scripts': [
            'RLPipes = rlpipes.cli:cli',
        ],
    }
)
