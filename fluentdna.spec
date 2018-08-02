# -*- mode: python -*-
from PyInstaller.compat import is_darwin
import gooey
gooey_root = os.path.dirname(gooey.__file__)
gooey_languages = Tree(os.path.join(gooey_root, 'languages'), prefix = 'gooey/languages')
gooey_images = Tree(os.path.join(gooey_root, 'images'), prefix = 'gooey/images')

block_cipher = None
pathextras=['D:\\josiah\\Projects\\DDV']
mainexepath=['DDV\\fluentdna.py']
excludelibs=[]
if is_darwin: pathextras=[] 
if is_darwin: excludelibs=['FixTk', 'tcl', 'tk', '_tkinter', 'tkinter', 'Tkinter']
if is_darwin: mainexepath=['DDV/fluentdna.py']

example_data=[('DDV/example_data', 'DDV/example_data'),]
a = Analysis(mainexepath,
   pathex=pathextras,
   binaries=[],
   datas=[('DDV/html_template','DDV/html_template')],
   hiddenimports=[],
   hookspath=[],
   runtime_hooks=[],
   excludes=excludelibs,
   win_no_prefer_redirects=False,
   win_private_assemblies=False,
   cipher=block_cipher)
             
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
options = [('u', None, 'OPTION'), ('u', None, 'OPTION'), ('u', None, 'OPTION')]

exe = EXE(pyz,
   a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          options,
          gooey_languages, # Add them in to collected files
          gooey_images, # Same here.

   name='FluentDNA',
   debug=False,
   strip=None, #strip=False,
   upx=False, #executable compression
   console=False,
   windowed=True,
   icon=os.path.join(gooey_root, 'images', 'program_icon.ico'))
   #(gooey_root, 'DDV','html_template','img','favicon.ico'))


#coll = COLLECT(exe,
#   a.binaries,
#   a.zipfiles,
#   a.datas,
#   strip=False,
#   upx=True,
#   name='FluentDNA')
   
#the app will be useful ONLY for the MAC GUI/windowed version
#if is_darwin:
#   app = BUNDLE(exe,
#                name='FluentDNA.app',
#                icon='favicon.icns')
