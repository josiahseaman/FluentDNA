# -*- mode: python -*-
from PyInstaller.compat import is_darwin

block_cipher = None
pathextras=['D:\\josiah\\Projects\\DDV']
mainexepath=['DDV\\fluentdna.py']
excludelibs=[]
if is_darwin: pathextras=[]
if is_darwin: excludelibs=['FixTk', 'tcl', 'tk', '_tkinter', 'tkinter', 'Tkinter']
if is_darwin: mainexepath=['DDV/fluentdna.py']

a = Analysis(mainexepath,
   pathex=pathextras,
   binaries=[],
   datas=[('DDV/example_data', 'DDV/example_data'),
   ('DDV/html_template','DDV/html_template')],
   hiddenimports=[],
   hookspath=[],
   runtime_hooks=[],
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
