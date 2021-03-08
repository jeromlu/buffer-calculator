# -*- mode: python -*-

block_cipher = None


a = Analysis(['buffer_calculator.pyw'],
             pathex=['C:\\Users\\jeromlu2\\LukaFiles\\04_Programiranje\\01_Python\\02_MyProjects\\buffer-calculations\\pH_prediction'],
             binaries=[],
             datas=[('./testing/calc_pH_values.exe', './testing/'),
					('./testing/version.dat', './testing/')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='buffer_calculator',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False,
		  icon ='c:\\Users\\jeromlu2\\LukaFiles\\04_Programiranje\\01_Python\\02_MyProjects\\buffer-calculations\\pH_prediction\\resources\\icons\\laboratory_con.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='buffer_calculator')
