rem This first line is necessary to run Cd1x_e.exe from the directory it is located in.
cd /D "%~dp0"
cls
del psplot.ps
del psplot.pdf
del psplot.log
Cd1x_e.exe
rem psplot.ps
