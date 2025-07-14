@echo off
setlocal enabledelayedexpansion

REM Create Feature Branch Script for Windows
REM Usage: scripts\create-feature.bat <component> <description>
REM Example: scripts\create-feature.bat equibind integration

REM Colors for output (Windows 10+)
set "BLUE=[94m"
set "GREEN=[92m"
set "YELLOW=[93m"
set "RED=[91m"
set "NC=[0m"

REM Function to print colored output
:print_status
echo %BLUE%[INFO]%NC% %~1
goto :eof

:print_success
echo %GREEN%[SUCCESS]%NC% %~1
goto :eof

:print_warning
echo %YELLOW%[WARNING]%NC% %~1
goto :eof

:print_error
echo %RED%[ERROR]%NC% %~1
goto :eof

REM Check if we're in a git repository
git rev-parse --git-dir >nul 2>&1
if errorlevel 1 (
    call :print_error "Not in a git repository. Please run this script from the repository root."
    exit /b 1
)

REM Check arguments
if "%~2"=="" (
    call :print_error "Usage: %0 <component> <description>"
    call :print_error "Example: %0 equibind integration"
    call :print_error "Example: %0 gnina enhancement"
    exit /b 1
)

set "COMPONENT=%~1"
set "DESCRIPTION=%~2"

REM Convert to lowercase and replace spaces with hyphens
for /f "tokens=*" %%i in ('echo %COMPONENT% ^| powershell -Command "$input = [Console]::In.ReadLine(); $input.ToLower() -replace ' ', '-'"') do set "COMPONENT=%%i"
for /f "tokens=*" %%i in ('echo %DESCRIPTION% ^| powershell -Command "$input = [Console]::In.ReadLine(); $input.ToLower() -replace ' ', '-'"') do set "DESCRIPTION=%%i"

set "BRANCH_NAME=feature/%COMPONENT%-%DESCRIPTION%"

call :print_status "Creating feature branch: %BRANCH_NAME%"

REM Check if branch already exists locally
git show-ref --verify --quiet refs/heads/%BRANCH_NAME%
if not errorlevel 1 (
    call :print_error "Branch %BRANCH_NAME% already exists locally"
    exit /b 1
)

REM Check if branch exists on remote
git ls-remote --heads origin %BRANCH_NAME% | findstr /c:"%BRANCH_NAME%" >nul
if not errorlevel 1 (
    call :print_error "Branch %BRANCH_NAME% already exists on remote"
    exit /b 1
)

REM Fetch latest changes
call :print_status "Fetching latest changes from remote..."
git fetch origin

REM Check if develop branch exists
git show-ref --verify --quiet refs/heads/develop
if errorlevel 1 (
    call :print_warning "Develop branch doesn't exist. Creating from main..."
    git checkout main
    git pull origin main
    git checkout -b develop
    git push -u origin develop
)

REM Switch to develop and update
call :print_status "Updating develop branch..."
git checkout develop
git pull origin develop

REM Create feature branch
call :print_status "Creating feature branch from develop..."
git checkout -b %BRANCH_NAME%

REM Push to remote
call :print_status "Pushing feature branch to remote..."
git push -u origin %BRANCH_NAME%

call :print_success "Feature branch '%BRANCH_NAME%' created successfully!"
call :print_status "You can now start developing your feature."
call :print_status "When ready, create a Pull Request from %BRANCH_NAME% to develop"

REM Show next steps
echo.
call :print_status "Next steps:"
echo   1. Start developing your feature
echo   2. Make regular commits with descriptive messages
echo   3. Push changes: git push origin %BRANCH_NAME%
echo   4. Create Pull Request to develop when ready
echo   5. Clean up branch after merge: git branch -d %BRANCH_NAME%

endlocal 