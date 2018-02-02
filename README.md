# FluentDNA Data Visualization tool (DDV) 

This application creates visualizations of FASTA formatted DNA nucleotide data.
FluentDNA generates a DeepZoomImage visualizations similar to Google Maps for FASTA files.

## Quick Start
You will need:
1. Python 3.4 or older: [Download Link](https://www.python.org/downloads/release/python-343/)
2. Git: [Download Link](https://git-scm.com/downloads)

**Installation**  
From a command line in your python directory or virtual environment:    
`pip install --process-dependency-links git+https://github.com/josiahseaman/FluentDNA.git@pip`

**Running**  
fluentdna.py will be placed in the scripts folder and accessible through PYTHONPATH.

`python /path/to/site-packages/DDV/fluentdna.py --fasta="/path/to/site-packages/DDV/example_data/hg38_chr19_sample.fa"`  
**Note:** Since Windows ignores the #!/bin/usr/python line, you'll need to use python and the full path to the script:  
`python C:\yourvenv\Scripts\fluentdna.py --fasta="C:\path\to\yourfasta.fa"`  

To use the interactive browser, especially for large files, start a server.  

`python /path/to/site-packages/DDV/fluentdna.py --runserver`  
Then open your browser and enter the URL: `http://localhost:8000/`  

To run FluentDNA from your own python script I recommend looking at fluentdna.py for examples such as `create_tile_layout_viz_from_fasta()`

### Example Commands
Check out the file [example_DDV_commands.txt](https://github.com/josiahseaman/DDV/blob/python-master/example_DDV_commands.txt) for more examples.

## Support Contact 
If you run into any problems or would like to use DDV in research, contact me at **josiah@newline.us**.  I'm happy to support my own software and always interested in new collaborations.

## DDV 2.0 Features

DDV 2.0 is a complete rewrite in Python of DDV.  DDV 2.0 has a much expanded feature set for handling
large, multipart files.  It can put an entire genome on a single image, marked with contig names.
DDV 2.0 has features in development for exploring genome alignments, annotations, and transposon alignments.
It was developed by Newline Technical Innovations and can be found at:
https://github.com/josiahseaman/DDV/tree/python-master


# Compile Instructions
Intended for developers:

### Windows:

  - Install Mercurial so it is properly on your path
  - Download PyWin32 from: http://sourceforge.net/projects/pywin32/files/pywin32/Build%20219/pywin32-219.win-amd64-py3.4.exe/download
  - Using the easy_install in your new Virtual Environment (/path/to/adsm_venv/Scripts/easy_install), install PyWin32:

        `easy_install pywin32-219.win-amd64-py3.4.exe`

  - Using /path/to/adsm_venv/Scripts/pip install cx_freeze:

        `pip install hg+https://bitbucket.org/BryanHurst/cx_freeze`

    > Note: this currently does not work; instead, install cx_freeze 4.3.4 using pip, then manually apply [this patch](https://bitbucket.org/BryanHurst/cx_freeze/commits/eba6cb644d390f69f07adbf9fdcead71ec0feebf?at=default) and [this patch](https://bitbucket.org/BryanHurst/cx_freeze/commits/22d73fe6386d92834339bdea30b3786a3543b2de?at=default) to the cx_freeze files that pip installed in your site-packages folder.

### Linux:

  - Requires ldd and objdump installed (probably already on your system)
  - Install Mercurial `sudo apt-get install mercurial`
  - You need a custom compiled version of Python3.4 (will use instead of venv)

        sudo apt-get install zlib1g-dev libbz2-dev libncurses5-dev libreadline6-dev libsqlite3-dev libssl-dev libgdbm-dev liblzma-dev tk8.5-dev
        wget https://www.python.org/ftp/python/3.4.3/Python-3.4.3.tgz
        tar zxvf Python-3.4.3.tgz
        rm Python-3.4.3.tgz
        cd Python-3.4.3/
        ./configure --prefix=/path/to/projects/ddv_python --exec_prefix=/path/to/projects/ddv_python
        make
        make altinstall
        /path/to/projects/ddv_python/bin/pip uninstall setuptools
        /path/to/projects/ddv_python/bin/pip uninstall pip
        wget https://pypi.python.org/packages/source/s/setuptools/setuptools-3.4.4.tar.gz
	    tar -vzxf setuptools-3.4.4.tar.gz
	    rm setuptools-3.4.4.tar.gz
	    cd setuptools-3.4.4
        /path/to/projects/ddv_python/bin/python setup.py install
        cd ..
        rm -r setuptools-3.4.4/
        wget https://pypi.python.org/packages/source/p/pip/pip-1.5.6.tar.gz
	    tar -vzxf pip-1.5.6.tar.gz
	    rm pip-1.5.6.tar.gz
	    cd pip-1.5.6
	    /path/to/projects/ddv_python/bin/python setup.py install
	    cd ..
	    rm -r pip-1.5.6

  - Using the new python, install all the requirements `/path/to/projects/ddv_python/bin/pip install -r /path/to/DDV/Requirements.txt`
  - `/path/to/projects/ddv_python/bin/pip install hg+https://bitbucket.org/BryanHurst/cx_freeze`
    - If the above install fails, then there is a problem with your python shared libraries, I have a clone of the cx_freeze repo with a temp fix
      - CD to a directory where you want to download it, then `hg clone hg+https://bitbucket.org/BryanHurst/cx_freeze; cd cx_freeze; /path/to/projects/adsm_python/bin/python setup.py install`

## History: DDV 1.1
This project is a fork of the C# DDV developed at Concordia University.
https://github.com/photomedia/DDV/

DDV Licence:
https://github.com/photomedia/DDV/blob/master/DDV-license.txt

### Examples (Demonstration):

http://www.photomedia.ca/DDV/

Visualizations generated with DDV can be placed on a web server. 
The following contains the links to examples of the visualizations 
generated by Tomasz Neugebauer, Éric Bordeleau, Vincent Burrus and Ryszard Brzezinski 
with this software: DNA Data Visualizations Generated with DDV Software. 
These examples include a number of bacteria chromosomes, as well as the entire Homo Sapiens genome. 
