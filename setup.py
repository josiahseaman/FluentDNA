from DDV.DDV import __version__

setup(name='DDV',
      version=__version__,
      description='DDV Application',
      options={'build_exe': build_exe_options,
               'install_exe': {'build_dir': build_exe_options['build_exe']}},
      executables=[Executable('DDV.py', base=base, icon='favicon.ico', targetName='DDV'+EXTENSION), ],
      cmdclass=cmdclass,
      install_requires=requirements,
      dependency_links=urls
      )


from setuptools import setup, find_packages

setup(
    name='DDV',
    version=__version__,
    description='Visualization tool for fasta files.  Supports whole genome alignment and multiple sequence alignment.',
    author='Josiah Seaman, Bryan Hurst',
    author_email='josiah.seaman@gmail.com',
    license='BSD',
    packages=['DDV'],  # TODO: revise this line
    include_package_data=True,
    install_requires=[
        'Pillow==3.2.0',
        'six==1.10.0',
        'wheel==0.24.0',
        'psutil==4.3.1',
        'blist==1.3.6',
        'numpy==1.13.3',
        'natsort==5.1.1',
    ],
    url='https://github.com/josiahseaman/DDV',
    download_url='https://github.com/josiahseaman/DDV',  # TODO: post a tarball
    keywords=['bioinformatics', 'dna', 'fasta', 'chain', 'alignment'],
    classifiers=[
        'Development Status :: 4 - Beta',  # 5 - Production/Stable
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)