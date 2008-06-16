for ext in txt html css png js doctree pickle json
do
    find build/ -name "*.${ext}" -exec rm -v {} \;
done

PYTHONPATH=$HOME/Development/yt/trunk/
make html

echo "Files deleted since last svn commit:"
echo 
svn stat | grep \!
echo

echo "Files added since last svn commit:"
echo " (generated with 'svn stat | grep \? )'"
echo
svn stat | grep \?
echo
