from setuptools import setup, find_packages

def readme():
    with open("README.md", 'r') as f:
    	 return f.read()

setup(
	name='symdr',
	version='1.0.0',
	author='A.Dzhanbekova, S.Kotov, M.Malyutin, M.Samoilov, V.Utupyina, M.Arendarenko, T.Savvateeva',
        author_email='',
	description='Find dispersion relation in PDE, systems of PDE, and discrete analogs',
	long_description=readme(),
	long_description_content_type='text/markdown',
	url='',
	packages=find_packages(),
  	install_requires=['sympy>=1.12'],
 	 classifiers=[
    	 'Programming Language :: Python :: 3.11',
    	 'License :: OSI Approved :: MIT License',
    	 'Operating System :: OS Independent'
  	 ],
 	  keywords='PDE, Continuum mechanics, computer algebra, python',
 	  project_urls={
    	 'github': 'https://github.com/symdr/symdr'
  	 },
 	  python_requires='>=3.7'

)
