@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul
cd /d "%~dp0"
echo FermatOrbits: проверка результатов статьи (орбиты K_{n,p} и полный поиск в GL_3)
echo Аргументы: [orbits NMAX PMAX ^| gl3 p [p ...]] [-1]; без аргументов - полный протокол
echo.

set "SDKMAJOR="
where dotnet >nul 2>&1
if not errorlevel 1 (
    for /f "tokens=1 delims=." %%v in ('dotnet --version 2^>nul') do set "SDKMAJOR=%%v"
)
echo !SDKMAJOR!| findstr /r "^[0-9][0-9]*$" >nul || set "SDKMAJOR="

if defined SDKMAJOR (
    echo Найден .NET SDK !SDKMAJOR!.x - сборка для net!SDKMAJOR!.0 ^(Release^)
    call :writeproj net!SDKMAJOR!.0
    if exist obj rd /s /q obj
    if exist bin rd /s /q bin
    if exist build rd /s /q build
    dotnet build FermatOrbits.csproj -c Release -nologo -v q -o build
    if errorlevel 1 goto :fail
    echo.
    build\FermatOrbits.exe %*
    if errorlevel 1 goto :fail
) else (
    echo .NET SDK не найден - используется компилятор csc.exe из .NET Framework 4.x
    set "CSC=%SystemRoot%\Microsoft.NET\Framework64\v4.0.30319\csc.exe"
    if not exist "!CSC!" set "CSC=%SystemRoot%\Microsoft.NET\Framework\v4.0.30319\csc.exe"
    if not exist "!CSC!" (
        echo Компилятор не найден. Установите .NET SDK: https://dotnet.microsoft.com/download
        goto :fail
    )
    "!CSC!" /nologo /optimize+ /define:LEGACY /out:FermatOrbits.exe Program.cs
    if errorlevel 1 goto :fail
    echo.
    FermatOrbits.exe %*
    if errorlevel 1 goto :fail
)

echo.
echo Протокол записан в файл FermatOrbits_protocol.txt в этой папке
goto :end

:fail
echo.
echo Сборка или запуск завершились с ошибкой.

:end
pause
exit /b

:writeproj
rem Файл проекта создаётся заново под установленный SDK: %~1
(
echo ^<Project Sdk="Microsoft.NET.Sdk"^>
echo   ^<PropertyGroup^>
echo     ^<OutputType^>Exe^</OutputType^>
echo     ^<TargetFramework^>%~1^</TargetFramework^>
echo     ^<RollForward^>LatestMajor^</RollForward^>
echo     ^<AssemblyName^>FermatOrbits^</AssemblyName^>
echo     ^<RootNamespace^>FermatOrbits^</RootNamespace^>
echo     ^<Nullable^>disable^</Nullable^>
echo     ^<ImplicitUsings^>disable^</ImplicitUsings^>
echo     ^<Optimize^>true^</Optimize^>
echo     ^<NoWarn^>CA1416^</NoWarn^>
echo     ^<SatelliteResourceLanguages^>en^</SatelliteResourceLanguages^>
echo   ^</PropertyGroup^>
echo ^</Project^>
) > FermatOrbits.csproj
exit /b 0
