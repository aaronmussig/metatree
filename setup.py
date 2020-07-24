import os
import re

from setuptools import setup


def read_meta():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'metatree/__init__.py')
    with open(path) as fh:
        hits = re.findall(r'__(\w+)__ ?= ?["\'](.+)["\']\n', fh.read())
    return {k: v for k, v in hits}


def readme():
    with open('README.md') as f:
        return f.read()


meta = read_meta()

setup(name=meta['title'],
      version=meta['version'],
      description=meta['description'],
      author=meta['author'],
      author_email=meta['author_email'],
      url=meta['url'],
      license=meta['license'],
      long_description=readme(),
      long_description_content_type='text/markdown',
      project_urls={
          "Bug Tracker": "https://github.com/aaronmussig/PhyloDM/issues",
          "Documentation": "https://github.com/aaronmussig/PhyloDM",
          "Source Code": "https://github.com/aaronmussig/PhyloDM",
      },
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='phylogenetic distance matrix symmetric',
      packages=['metatree', 'metatree.io', 'metatree.external'],
      entry_points={
          'console_scripts': [
              'metatree = metatree.__main__:main'
          ]
      },
      install_requires=['phylorank>0.1.0', 'genometreetk>0.1.2', 'dendropy>=4.1.0', 'tqdm>=4.31.0', 'biolib>=0.1.0',
                        'biopython', 'seaborn', 'matplotlib', 'numpy', 'ete3', 'pandas', 'scipy'],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])]
      )
