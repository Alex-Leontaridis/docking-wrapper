@echo off
REM Dummy DiffDock script for testing on Windows

echo [DUMMY DIFFDOCK] Called with args: %*

REM Find output directory argument (after --out_dir)
set "out_dir="
set "prev=0"
for %%i in (%*) do (
    if "!prev!"=="1" (
        set "out_dir=%%i"
        goto :found_out
    )
    if "%%i"=="--out_dir" (
        set "prev=1"
    ) else (
        set "prev=0"
    )
)

:found_out
REM Create dummy output files if specified
if defined out_dir (
    if not exist "%out_dir%" mkdir "%out_dir%"
    echo REMARK DIFFDOCK OUTPUT > "%out_dir%\diffdock_out.sdf"
    echo REMARK DIFFDOCK CONFIDENCE > "%out_dir%\diffdock_confidence.txt"
)

exit /b 0 