@echo off
setlocal enabledelayedexpansion
REM Dummy Vina script for testing on Windows

echo [DUMMY VINA] Called with args: %*

REM Find output file argument (after --out)
set "out_file="
set "prev=0"
for %%i in (%*) do (
    if "!prev!"=="1" (
        set "out_file=%%i"
        goto :found_out
    )
    if "%%i"=="--out" (
        set "prev=1"
    ) else (
        set "prev=0"
    )
)

:found_out
REM Create dummy output file if specified
if defined out_file (
    for %%f in ("%out_file%") do set "out_dir=%%~dpf"
    if not exist "!out_dir!" mkdir "!out_dir!"
    echo REMARK VINA RESULT:    -4.145      0.000      0.000 > "%out_file%"
    echo MODEL        1 >> "%out_file%"
    echo ENDMDL >> "%out_file%"
)

exit /b 0 