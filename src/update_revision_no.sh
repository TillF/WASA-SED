#!/bin/sh
#update file General/svn_rev.f90 with current revision numbering to be included in the binary
#to be executed within the directory where you downloaded WASA-SED files into (manually or by Makefile)

# set output of git to variable
gitstr=`git describe --tags --long`

# check for local modifications
mod_files=`git status --porcelain`

if [ -z "$mod_files" ]
then
cmmt=""
else
cmmt="(local modifications)"
fi

REP_DATE=`git show -s --format=%ci`
DATE=` date +%Y-%m-%d" "%H:%M`
echo "!this file is updated on calling update_revision_no.sh/.bat to create svn_rev.txt which is linked to the source code
rev_string1='$gitstr $cmmt'
rev_string2='repository date $REP_DATE, built $DATE'" > $1/General/svn_rev.f90