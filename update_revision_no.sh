#!/bin/sh
#update file General/svn_rev.var with current revision numbering to be included in the binary
REV_NO=`svnversion`
if [ -z `echo $REV_NO | grep M` ]
then
MODS=""
else
MODS="(mixed revisions , local modifications)"
fi

REP_DATE=`svn info | grep --only-matching "[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9][^(]*"`
DATE=`date`
echo "!this file is updated on calling update_revision_no.sh to create svn_rev.txt which is linked to the source code
rev_string1='$REV_NO $MODS'
rev_string2='repository date $REP_DATE, built 2011/03/18 10:28:53'" > General/svn_rev.var