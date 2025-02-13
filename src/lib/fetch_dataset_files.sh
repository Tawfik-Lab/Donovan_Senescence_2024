#!/bin/bash

TARGET_DIR="../../build/datasets"
mkdir -p "$TARGET_DIR"

# List of GEO dataset accession numbers
DATASETS=("GSE154659" "GSE155622" "GSE249746")

# Download GEO datasets from NCBI FTP
for dataset in "${DATASETS[@]}"
do
    echo "Downloading files for ${dataset} dataset..."

    # Create a directory for each dataset within the target directory
    mkdir -p "$TARGET_DIR/$dataset/downloads"
    PREFIX=${dataset:0:6}
    
    FTP_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}nnn/${dataset}/suppl/"
    
    # Download all files from the supplementary folder into the respective dataset directory
    echo "Downloading files for ${dataset} from ${FTP_URL}..."
    wget -r -nH --cut-dirs=6 -np -A "*.tar.gz,*.gz,*.txt,*.csv" -P "$TARGET_DIR/$dataset/downloads" "$FTP_URL"
done

# North et al. 2019: Download supplementary data linked to the manuscript
echo "Downloading files for North_2019 dataset..."
mkdir -p "$TARGET_DIR/North_2019/downloads"

# Download and place in the target directory
wget -O "$TARGET_DIR/North_2019/downloads/North_2019_supplementary_data.zip" "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/brain/142/5/10.1093_brain_awz063/2/awz063_supplementary_data.zip?Expires=1733830096&Signature=5LeRcwi~GCx-ZcDnNUyr733HxTlN-LeaGdXYaDLOvMorkGY7S0jDChFBaxFVpjqJu3oxI4tpYaAAxw-EJ5ZBwanxtWC~5VnHHnyhjLxEh7lxDyaaL9mtgmxllUPunKrn98LOripZ42k4K91U5l2r8PTV~Yq3ELj2sq1k9rIKu9ZTmqCsCaIYeWKdMfAHjJuGtvZqYuG2A1QjorgbVwzGiFR6P-YxHF10wG~C2q8RMPEovo37wwNBAkkpdvFV5HOOzf1XpqYQkzvuJ63LRgF1iw~9-9QYTBkyWTLKDmmrBUqMHudv4BocqBn2XPYh~JKzCP7-yb6GfwFWlY2WjlYh4w__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA"

# Unzip the supplementary data into the respective directory
echo "Unzipping North_2019/downloads/North_2019_supplementary_data.zip"
unzip "$TARGET_DIR/North_2019/downloads/North_2019_supplementary_data.zip" -d "$TARGET_DIR/North_2019/downloads"
