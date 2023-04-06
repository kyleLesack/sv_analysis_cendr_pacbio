#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for i in {1..10}
do
   sort set_${i}.txt > set_${i}_sorted.txt; mv set_${i}_sorted.txt set_${i}.txt
done


