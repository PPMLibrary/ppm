#!/usr/bin/env bash

OK="\033[0;32mok\033[0m"
FAIL="\033[0;31mfail!\033[0m"
WARN="\033[0;33mwarnings\033[0m"

echo $1 >> $2

cmd="$1 2>&1"

output=`eval $cmd`

if [ $? -eq 0 ]
then
    if [ -z "$output" ]
    then
	# all ok
	printf "$OK\n"
    else
	# warning
	printf "$WARN\n"
	echo "$output" >> $2
    fi
else
    # error
    printf "$FAIL\n"
    echo "$output" >> $2
    printf "\n\033[0;31m%s\033[0m (check %s for details)\n\n%s\n\n" "$3" "$2" "$output"
    exit 1
fi
