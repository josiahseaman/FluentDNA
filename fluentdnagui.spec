from PyInstaller.compat import is_darwin

block_cipher = None
pathextras=['D:\\josiah\\Projects\\DDV']
mainexepath=['DDV\\fluentdnagui.py']
excludelibs=[]
if is_darwin: pathextras=[] 
if is_darwin: excludelibs=['FixTk', 'tcl', 'tk', '_tkinter', 'tkinter', 'Tkinter']
if is_darwin: mainexepath=['DDV/fluentdnagui.py']

import gooey
gooey_root = os.path.dirname(gooey.__file__)
gooey_languages = Tree(os.path.join(gooey_root, 'languages'), prefix = 'gooey/languages')
gooey_images = Tree(os.path.join(gooey_root, 'images'), prefix = 'gooey/images')

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

options = [('u', None, 'OPTION')]

exe = EXE(pyz,
   a.scripts,
   exclude_binaries=True,
   name='FluentDNAgui',
   debug=False,
   strip=False,
   upx=True,
   console=False,
   icon='favicon.ico')

coll = COLLECT(exe,
   a.binaries,
   a.zipfiles,
   a.datas,
   gooey_languages, 
   gooey_images,
   strip=False,
   upx=True,
   name='FluentDNAgui')
          
          
#the app will be useful ONLY for the MAC GUI/windowed version
if is_darwin:
   app = BUNDLE(exe,
   				coll,
                name='FluentDNAgui.app',
                icon='favicon.icns')
