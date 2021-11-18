#!/bin/bash
file=para.txt
echo $1
n=0
while IFS= read -r line
do
     if test "$n" == "$1"
     then
            ./decomp $1 $line
     fi
     n=$((n+1))
done < "$file"
