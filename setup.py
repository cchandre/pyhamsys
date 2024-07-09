import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyhamsys',
    version='0.32',
    description='Some tools for Hamiltonian systems',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    url='http://github.com/cchandre/pyhamsys',
    classifiers=[
      'Programming Language :: Python :: 3',
      'Intended Audience :: Education',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: BSD License',
      'Operating System :: OS Independent',
      'Topic :: Scientific/Engineering',
      'Topic :: Scientific/Engineering :: Mathematics',
      'Topic :: Scientific/Engineering :: Chemistry',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Scientific/Engineering :: Astronomy'
    ],
    python_requires='>=3.6',
    author='Cristel Chandre',
    author_email='cristel.chandre@cnrs.fr',
    license='BSD',
    py_modules=['pyhamsys'],
    package_dir={'': 'pyhamsys/src'}, 
    install_requires=[
      "numpy",
      "scipy"]
)