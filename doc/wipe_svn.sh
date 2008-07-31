export PYTHONPATH=`pwd`/../
make html

echo "Files deleted since last svn commit:"
echo " (generated with 'svn stat | grep \! ')"
echo 
svn stat | grep \! |grep -v build
echo

echo "Files added since last svn commit:"
echo " (generated with 'svn stat | grep \? ')"
echo
svn stat | grep \? | grep -v build/html
echo

echo "Examples files status:"
echo
svn stat ../examples/

zip -r docs_html.zip build/html/
