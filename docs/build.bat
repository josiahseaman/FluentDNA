PyInstaller fluentdna.spec --onedir --runtime-hook=".\DDV\use_lib.py" --clean

:: you'll need to pip install pyinstaller==3.3.1
:: see installation.md for details
:: simpler: PyInstaller fluentdna.spec --clean --noconfirm