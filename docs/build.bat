PyInstaller fluentdna.spec --clean --noconfirm --onedir --runtime-hook=".\DDV\use_lib.py"

:: you'll need to pip install pyinstaller==3.3.1
:: see installation.md for details
:: simpler: PyInstaller fluentdna.spec --clean --noconfirm