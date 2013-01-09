#!/bin/sh

# OptionMatrix website Publishing.
# Copyright (C) 2012 Anthony Bradford.
echo
echo "Copyright (C) 2012 Anthony Bradford."

rm -f ftp.error

if [ "$#" != 3 ]; then
    echo
    echo
    echo "From the distribution top level directory run:"
    echo
    echo "./configure ; make website"
    echo
    echo
    exit 1
fi
echo
echo
echo Push doc/html to website
echo Push doc/html/optionmatrix.zip to website
echo Push $2 to website
echo Push doc/$3 to website
echo
echo The HTML contains hard links to http://opensourcefinancialmodels.com
echo
echo If you want to host the .pdf, tar.gz and .zip files you have to change
echo the hard links in doc/optionmatrix.texi.
echo
echo See:
echo "grep "http://opensourcefinancial" doc/optionmatrix.texi | egrep '.pdf|.zip|.gz'"
echo

#HOST='opensourcefinancialmodels.com'
echo

read -p "FTP website  : " HOST
read -p "FTP login    : " USER
read -p "FTP password : " PASSWD
echo

PUBLISH_TAR_GZ=0
PUBLISH_HTML=0
PUBLISH_PDF=0

while true; do 
    read -p "Publish OptionMatrix HTML documentation to the web (y/n)? " yn 
    case $yn in 
        [Yy]* ) PUBLISH_HTML=1 break;; 
        [Nn]* ) break;; 
        * ) echo "Please answer yes or no.";; 
    esac 
done 
if [ "$PUBLISH_HTML" = "1" ]; then
    make html

echo
echo
echo This might take some time. 
echo Directory creation errors are normal after initial run.
echo
echo

ftp -n $HOST <<END_SCRIPT 2>> ftp.error
quote USER $USER
quote PASS $PASSWD
binary
cd public_html
mkdir om
cd om
lcd doc/html
prompt
mput *.html
mkdir doc
cd doc
mput *.zip
cd ..
mkdir images
cd images
lcd images
mput *.png
quit
END_SCRIPT

fi

while true; do 
    echo
    read -p "Publish optionmatrix.tar.gz to the web (y/n)? " yn
    case $yn in 
        [Yy]* ) PUBLISH_TAR_GZ=1 break;; 
        [Nn]* ) break;; 
        * ) echo "Please answer yes or no.";; 
    esac 
done
if [ $PUBLISH_TAR_GZ = "1" ]; then
    make clean
    make dist

echo
echo
echo This might take some time.
echo Directory creation errors are normal after initial run.
echo
echo

ftp -n $HOST <<END_SCRIPT 2>> ftp.error
quote USER $USER
quote PASS $PASSWD
binary
cd public_html
mkdir om
cd om
mkdir downloads
cd downloads
put $2
mkdir rpm
mkdir debian
cd debian
mkdir amd64
mkdir i386
quit
END_SCRIPT

fi

while true; do 
    echo
    read -p "Publish optionmatrix.pdf the web (y/n)? " yn
    case $yn in 
        [Yy]* ) PUBLISH_PDF=1 break;; 
        [Nn]* ) break;; 
        * ) echo "Please answer yes or no.";; 
    esac 
done 
if [ "$PUBLISH_PDF" = "1" ]; then
    make pdf

echo
echo
echo This might take some time.
echo Directory creation errors are normal after initial run.
echo
echo

ftp -n $HOST <<END_SCRIPT 2>> ftp.error
quote USER $USER
quote PASS $PASSWD
binary
cd public_html
mkdir om
cd om
mkdir doc
cd doc
lcd doc
put $3
quit
END_SCRIPT

fi

echo done

if [ -e ftp.error ]; then
    cat ftp.error
    rm -r ftp.error
fi
