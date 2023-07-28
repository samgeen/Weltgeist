!#/bin/bash
cd StellarSources
export FILE=StellarSources.zip
if [ -f "$FILE" ]; then
    echo "$FILE exists, to download again delete it first."
else 
    echo "$FILE does not exist, downloading."
    wget -O $FILE https://zenodo.org/record/8192816/files/StellarSources.zip?download=1
fi
unzip $FILE

