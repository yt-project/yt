find . -name "*.so" -exec rm -v {} \;
find . -name "*.pyc" -exec rm -v {} \;
find . -name "__config__.py" -exec rm -v {} \;
rm -rvf build dist
