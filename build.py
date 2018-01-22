import sys
import os
import stat
import pip
import shutil
import subprocess

from cx_Freeze import setup, Executable, build_exe
from importlib import import_module

from DDV import VERSION


BASE_DIR = os.path.dirname(os.path.abspath(__file__))

if sys.platform == 'win32':
    OS_DIR = 'windows'
    EXTENSION = '.exe'
    SCRIPT = '.cmd'
else:
    OS_DIR = 'linux'
    EXTENSION = ''
    SCRIPT = ''


def is_exe(file_path):
    access_mode = os.F_OK | os.X_OK
    if os.path.isfile(file_path) and not file_path.endswith('.bat') and not file_path.endswith('.sh') and os.access(file_path, access_mode):
        filemode = os.stat(file_path).st_mode
        ret = bool(filemode & stat.S_IXUSR or filemode & stat.S_IXGRP or filemode & stat.S_IXOTH)
        return ret


build_exe_options = {
    'build_exe': 'build',
    'optimize': 2,
    'excludes': [
        'PyInstaller',
        'tkinter',
        # CHANGE ME for any python packages in your project that you want excluded
    ],
    'includes': [
        # CHANGE ME for any python packages in your project that you want included
    ],
    'packages': [
        # Known missing Python imports.
        # This may need to be updated from time to time as new projects uncover more missing packages.
        'html',
        'shutil',
    ],
    'replace_paths': [('*', '')],
    'include_files': [
        # Standard Django items to bring in
        ('npu.exe' if sys.platform == 'win32' else 'npu', 'npu.exe' if sys.platform == 'win32' else 'npu'),  # Newline's Updater application for Django Desktop Core
        ('README.md', 'README.md'),

        # CHANGE ME for any files/folders you want included with your project
        ('DDV/html_template'),
    ],
    'include_msvcr': True  # CHANGE ME depending on if your project has licensing that is compatible with Microsoft's redistributable license
}
files = (file for file in os.listdir(BASE_DIR) if os.path.isfile(os.path.join(BASE_DIR, file)))
for file in files:
    # We need to check for the npu here instead of just including it as it would come up as an exe in the is_exe check
    if [file for part in ['.so', '.dll', '.url', 'webpack-stats.json'] if part.lower().split(' ')[0] in file.lower()] or is_exe(os.path.join(BASE_DIR, file)) and 'npu' not in file:
        build_exe_options['include_files'].append((file, file))


def query_yes_no(question, default='yes'):
    valid = {'yes': True, 'y': True, "no": False, 'n': False}

    if default is None:
        prompt = " [y/n] "
    elif default in ['yes', 'y']:
        prompt = " [Y/n] "
    elif default in ['no', 'n']:
        prompt = " [y/N] "
    else:
        raise ValueError("Invalid default answer!")

    while True:
        sys.stdout.write('\n' + question + prompt)

        choice = input().lower()

        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no'.\n")


def get_packages_in_path(path):
    packages = []
    for folder in os.listdir(path):
        if os.path.exists(os.path.join(path, folder, '__init__.py')):
            packages.extend([folder, ])
        elif folder.endswith('.egg') and not os.path.isfile(os.path.join(path, folder)):
            for subfolder in os.listdir(os.path.join(path, folder)):
                if os.path.exists(os.path.join(path, folder, subfolder, '__init__.py')):
                    packages.extend([subfolder, ])
    return packages


def get_installed_packages():
    site_packages = pip.util.site_packages if hasattr(pip, 'util') else pip.utils.site_packages

    return get_packages_in_path(site_packages)


def get_local_packages():
    return get_packages_in_path(BASE_DIR)


def remove_empty_folders(path):
    if not os.path.isdir(path):
        return

    files = os.listdir(path)
    if len(files):
        for f in files:
            fullpath = os.path.join(path, f)
            if os.path.isdir(fullpath):
                remove_empty_folders(fullpath)

    files = os.listdir(path)
    if len(files) == 0:
        print("Removing empty folder:", path)
        os.rmdir(path)


def parse_requirements_and_links(requirements_file, existing_requirements=None, existing_links=None):
    """
    Proper Git lines from Requirements.txt:
        git+https://git.myproject.org/MyProject.git
        git+https://git.myproject.org/MyProject.git@v1.0.1

    Proper Mercurial lines from Requirements.txt:
        hg+https://hg.myproject.org/MyProject/
        hg+https://hg.myproject.org/MyProject/#egg=MyProject
        hg+https://hg.myproject.org/MyProject/@v1.0.1#egg=MyProject
    """
    if not existing_requirements:
        existing_requirements = []
    if not existing_links:
        existing_links = []

    with open(requirements_file, 'r') as requirements:
        for line in requirements:
            line = line.strip()

            version = None
            package = None
            link = None

            if line.startswith('git+'):
                parts = line.split('@')
                if parts.__len__() > 1:
                    version = parts[1]

                url_parts = parts[0].split('/')
                package = url_parts[-1].split('.git')[0] if '.git' in url_parts[-1] else None

                if version:
                    link = line + "#egg=" + package + "-" + version
                    if package:
                        package = package + '==' + version
                else:
                    link = line
            elif line.startswith('hg+'):
                line = line.split('#egg=')[0]

                parts = line.split('@')
                url_parts = parts[0].split('/')
                package = url_parts[-1]
                if parts.__len__() > 1:
                    version = parts[1]

                if version:
                    link = line + '#egg=' + package + '-' + version
                    if package:
                        package = package + '==' + version
                else:
                    link = line
            else:
                if line:
                    package = line

            if package:
                existing_requirements.append(package)
                if link:
                    existing_links.append(link)

    return existing_requirements, existing_links


# TODO: This build doesn't quite work in Linux. The files don't end up in the proper location.
# Currently you must manually move files around after the build finishes to get it to work in Linux.
class BuildDDV(build_exe):
    def run(self):
        print("\nYou should only run this build script if you are a CLEAN VirtualEnv!\n"
              "The VirtualEnv will be deployed with the project, "
              "so make sure it ONLY has the REQUIRED dependencies installed!")
        if not query_yes_no("Are you in a CLEAN Python Environment?", default='no'):
            sys.exit()

        shutil.rmtree(os.path.join(BASE_DIR, self.build_exe), ignore_errors=True)

        # Grab all the packages that we should include (local and those installed in the virtualenv)
        self.packages.extend(get_installed_packages())
        self.packages.extend(get_local_packages())
        # Cleanup packages to account for any excludes
        self.packages = [package for package in self.packages if package not in self.excludes]

        # Check all installed packages and see if they list any dependencies they need copied when frozen
        for package in self.packages:
            try:
                package = import_module(package)
                if getattr(package, '__include_files__', False):
                    self.include_files.extend(package.__include_files__)
            except Exception as e:
                print("Error bringing in dependent files!", str(e), "This is probably okay.")
                continue

        build_exe.run(self)

        # Cleanup the build dir
        if sys.platform != 'win32':  # If we are on a unix type system
            lib_files = os.listdir(os.path.join(BASE_DIR, self.build_exe, 'lib'))
            for file in lib_files:
                if "python" not in str(file).lower():
                    shutil.move(os.path.join(BASE_DIR, self.build_exe, 'lib', file), os.path.join(BASE_DIR, self.build_exe, file))
            if os.path.exists(os.path.join(BASE_DIR, self.build_exe, 'bin', 'env')):
                env_files = os.listdir(os.path.join(BASE_DIR, self.build_exe, 'bin', 'env'))
                for file in env_files:
                    shutil.move(os.path.join(BASE_DIR, self.build_exe, 'bin', 'env', file), os.path.join(BASE_DIR, self.build_exe, file))  

        # move binary dependencies into bin/env
        files = (file for file in os.listdir(os.path.join(BASE_DIR, self.build_exe)) if os.path.isfile(os.path.join(BASE_DIR, self.build_exe, file)))
        os.makedirs(os.path.join(BASE_DIR, self.build_exe, 'bin', 'env'))
        for file in files:
            # TODO: Check for linux python.so files
            if not [file for part in ['library.zip', 'README.md', 'python34.dll', 'MSVCR100.dll', 'npu', '.url', 'webpack-stats.json'] if part.lower().split(' ')[0] in file.lower()] and not is_exe(os.path.join(BASE_DIR, self.build_exe, file)):  #NOTE: The split here could cause issues and is speculative
                shutil.move(os.path.join(BASE_DIR, self.build_exe, file),
                            os.path.join(BASE_DIR, self.build_exe, 'bin', 'env', file))

        # Look for any DLLs that the included packages may have and copy them into bin/env as well
        for package in self.packages:
            try:
                package = import_module(package)
                location = os.path.dirname(package.__file__)
                for root, dirnames, filenames in os.walk(location):
                    for filename in filenames:
                        if filename.lower().split('.')[-1] in 'dll so'.split():
                            shutil.copy(os.path.join(root, filename), os.path.join(BASE_DIR, self.build_exe, 'bin', 'env', filename))
            except:
                continue


base = 'Console'
requirements, urls = parse_requirements_and_links(os.path.join(BASE_DIR, 'Requirements.txt'))

cmdclass = {"build_exe": BuildDDV,}

setup(name='FluentDNA',
      version=VERSION,
      description='DNA Sequence Visualization',
      options={'build_exe': build_exe_options,
               'install_exe': {'build_dir': build_exe_options['build_exe']}},
      executables=[Executable('fluentdna.py', base=base, icon='favicon.ico', targetName='fluentdna'+EXTENSION), ],
      cmdclass=cmdclass,
      install_requires=requirements,
      dependency_links=urls
      )

# Cleanup step after any sort of setup operation
# TODO: See if this causes issues at the end of an 'install' command.
# If so, override the msi build command and put this at the end of the normal build and the installer build
remove_empty_folders(os.path.join(BASE_DIR, build_exe_options['build_exe']))
