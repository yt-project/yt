for fn in *.py
do
    COUNT=`cat *.rst | grep --count ${fn}`
    if [ $COUNT -lt 1 ]
    then
        echo ${fn}
    fi
done
