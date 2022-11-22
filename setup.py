from setuptools import setup, find_packages
from pkg_resources import parse_requirements
from pathlib import Path

with Path('requirements.txt').open() as requirements_txt:
    install_requires = [
        str(requirement)
        for requirement
        in parse_requirements(requirements_txt)
    ]

# TODO: fix tests incorrectly running when installing with setup.py
setup(
    name="make_prg",
    version="0.4.0",
    packages=find_packages(),
    url="https://github.com/rmcolq/make_prg",
    license="MIT",
    entry_points={"console_scripts": ["make_prg = make_prg.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3", "hypothesis >= 4.0", "pytest"],
    install_requires=install_requires,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
