IF NOT DEFINED CIL_VERSION (
ECHO CIL_VERSION Not Defined.
exit 1
)

xcopy /e "%RECIPE_DIR%\..\..\..\*" "%SRC_DIR%"
cd %SRC_DIR%\Wrappers\Python
%PYTHON% setup.py install --single-version-externally-managed --record=record.txt
if errorlevel 1 exit 1
