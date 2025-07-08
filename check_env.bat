@echo off
echo Checking conda environment...

REM Add conda to PATH
set PATH=%PATH%;%USERPROFILE%\miniconda3;%USERPROFILE%\miniconda3\Scripts

echo Listing conda environments:
conda env list

echo.
echo Checking if equibind environment exists:
conda env list | findstr equibind

echo.
echo Listing packages in equibind environment:
conda list -n equibind

pause 