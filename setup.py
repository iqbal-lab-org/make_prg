from setuptools import setup, find_packages

setup(
    name="make_prg",
    version="0.0.0",
    packages=find_packages(),
    url="https://github.com/rmcolq/make_prg",
    license="MIT",
    entry_points={"console_scripts": ["make_prg = make_prg.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=[
        "biopython>=1.70",
        "numpy>=1.14.0",
        "scikit-learn>=0.19.1",
        "scipy>=1.0.1",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
