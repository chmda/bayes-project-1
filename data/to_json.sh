#!/bin/bash
SCRIPT=convert.R

echo "Converting data to JSON file...";
Rscript $SCRIPT
echo "Done.";