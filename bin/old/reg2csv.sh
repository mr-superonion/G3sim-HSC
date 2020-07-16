#!/usr/bin/sh
sed '/dip/d;s/^circle(//;s/d) #.*$//' $1
