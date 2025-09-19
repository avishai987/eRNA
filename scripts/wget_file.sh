#!/bin/bash
# Downloads a file from a given URL using wget.
# Resumes the download if interrupted ( -c).
# Retries indefinitely until the download is complete(--tries=0).
#good for large files like fasta files.
# Usage: ./wget_file.sh <URL> <DIRECTORY>

url=$1
dir=$2
cd $dir
wget -c --tries=0 $url