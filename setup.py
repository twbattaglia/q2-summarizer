import ast
from setuptools import setup, find_packages
import re

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2_summarizer/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

setup(
    name="q2-summarizer",
    version=version,
    packages=find_packages(),
    install_requires=['qiime >= 2.0.0',
                      'pandas', 'biom-format',
                      'seaborn',
                      'jinja2'],
    author="Thomas W. Battaglia",
    author_email="tb1280@nyu.edu",
    description="QIIME2 plugin for generating interactive summary data",
    entry_points={
        "qiime.plugins":
        ["q2-summarizer=q2_summarizer.plugin_setup:plugin"]
    }
)
