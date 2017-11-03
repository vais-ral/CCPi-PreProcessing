if [ -z "$CIL_VERSION" ]; then
    echo "Need to set CIL_VERSION"
    exit 1
fi  
cp -r "${RECIPE_DIR}/../../" ${SRC_DIR}
cd ${SRC_DIR}/Python
$PYTHON setup.py install
