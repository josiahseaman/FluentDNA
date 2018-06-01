# -*- mode: python -*-

block_cipher = None

# Icons: https://stackoverflow.com/questions/9946760/add-image-to-spec-file-in-pyinstaller
a = Analysis(['DDV\\fluentdna.py'],
             pathex=['D:\\josiah\\Projects\\DDV'],
             binaries=[],
             datas=[('DDV/example_data', 'DDV/example_data'),
             ('DDV/html_template','DDV/html_template')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
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
          strip=False,
          upx=True,
          console=True,
          icon='D:\\josiah\\Projects\\DDV\\favicon.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='FluentDNA')