# -*- mode: python -*-
from PyInstaller.compat import is_darwin
import os, sys

block_cipher = None
pathextras = os.path.dirname(sys.argv[0])
mainexepath = [os.path.join('FluentDNA','fluentdna.py')]
excludelibs = []
if is_darwin:
    excludelibs = ['FixTk', 'tcl', 'tk', '_tkinter', 'tkinter', 'Tkinter']
    pathextras = os.path.dirname(pathextras)

a = Analysis(mainexepath,
   pathex=[pathextras],
   binaries=[],
   datas=[('FluentDNA/example_data', 'example_data'),
   ('FluentDNA/html_template','html_template'),
   ('docs','docs')],
   hiddenimports=['xml', 'pyexpat'],
   hookspath=[],
   runtime_hooks=['FluentDNA/use_lib.py'],
   excludes=excludelibs,
   win_no_prefer_redirects=False,
   win_private_assemblies=False,
   cipher=block_cipher)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
   a.scripts,
   exclude_binaries=True,
   name='FluentDNA',
   debug=False,
   strip=None, #strip=False,
   upx=True,
   console=True,
   icon='favicon.ico')

coll = COLLECT(exe,
   a.binaries,
   a.zipfiles,
   a.datas,
   strip=False,
   upx=True,
   name='FluentDNA')

#the app will be useful ONLY for the MAC GUI/windowed version
#if is_darwin:
#   app = BUNDLE(exe,
#                name='FluentDNA.app',
#                icon='favicon.icns')
