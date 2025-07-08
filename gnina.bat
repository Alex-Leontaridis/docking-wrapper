@echo off
setlocal enabledelayedexpansion
REM Dummy GNINA script for testing on Windows

REM Check for --version flag
set version_flag=0
for %%i in (%*) do (
    if "%%i"=="--version" set version_flag=1
)
if !version_flag! == 1 (
    echo GNINA Dummy Version 1.0.0
    exit /b 0
)

echo [DUMMY GNINA] Called with args: %*

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
    
    REM Create realistic GNINA output
    echo REMARK DUMMY GNINA OUTPUT > "%out_file%"
    echo REMARK CNN SCORE: -8.123 > "%out_file%"
    echo REMARK CNN SCORE: -7.456 > "%out_file%"
    echo REMARK CNN SCORE: -6.789 > "%out_file%"
    echo MODEL        1 >> "%out_file%"
    echo ATOM      1  C   UNL     1      -0.123   1.456   2.789  0.00  0.00           C >> "%out_file%"
    echo ATOM      2  O   UNL     1       1.234  -0.567   3.890  0.00  0.00           O >> "%out_file%"
    echo ATOM      3  N   UNL     1      -1.345   0.678   1.234  0.00  0.00           N >> "%out_file%"
    echo ENDMDL >> "%out_file%"
    echo MODEL        2 >> "%out_file%"
    echo ATOM      1  C   UNL     1       0.456  -1.234   2.567  0.00  0.00           C >> "%out_file%"
    echo ATOM      2  O   UNL     1       2.123   0.789   1.456  0.00  0.00           O >> "%out_file%"
    echo ATOM      3  N   UNL     1      -2.567   1.890   0.123  0.00  0.00           N >> "%out_file%"
    echo ENDMDL >> "%out_file%"
    echo MODEL        3 >> "%out_file%"
    echo ATOM      1  C   UNL     1       1.789   0.123  -1.456  0.00  0.00           C >> "%out_file%"
    echo ATOM      2  O   UNL     1       3.456   2.123   0.789  0.00  0.00           O >> "%out_file%"
    echo ATOM      3  N   UNL     1      -3.789   2.456  -0.567  0.00  0.00           N >> "%out_file%"
    echo ENDMDL >> "%out_file%"
)

echo [DUMMY GNINA] Completed successfully
exit /b 0 