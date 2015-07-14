#!/bin/bash

#
# run using eg:
# watch --interval 5 monitor_system.sh
#

echo =====DF
df -h
echo

echo =====MDSTAT
cat /proc/mdstat
echo

echo =====W
w
echo

#echo =====WHO
#who
#echo

echo ====TOP
top -b -n 1 | head -n 30
echo
