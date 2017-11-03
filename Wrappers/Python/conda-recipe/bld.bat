IF NOT DEFINED CIL_VERSION (
ECHO CIL_VERSION Not Defined.
exit 1
)

xcopy /e "%RECIPE_DIR%\..\..\..\" "%SRC_DIR%"
cd %SRC_DIR%\Python
%PYTHON% setup.py install
if errorlevel 1 exit 1
