#!/bin/bash

# Retrieve file
file=$1

# Count lines 
awk 'END{print NR}' $file

