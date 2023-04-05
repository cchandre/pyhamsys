import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='pyhamsys',
      version='0.0.1',
      description='Some tools for Hamiltonian systems',
      long_description=long_description,
      long_description_content_type="text/markdown",
      packages=setuptools.find_packages(),
      url='http://github.com/cchandre/pyhamsys',
      python_requires='>=3.8',
      author='Cristel Chandre',
      author_email='cristel.chandre@cnrs.fr',
      license='BSD',
      ppy_modules=['pyhamsys'],
      package_dir={'':'pyhamsys/src'}, 
      install_requires=[]
      )