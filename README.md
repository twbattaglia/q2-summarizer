# q2-summarizer
[![Build Status](https://travis-ci.org/twbattaglia/q2-summarizer.svg?branch=master)](https://travis-ci.org/twbattaglia/q2-summarizer)  
QIIME2 plugin for generating interactive summary data  

## Installation
See https://github.com/qiime2/qiime2/wiki/Installing-and-using-QIIME-2 for more information on how to set up a qiime2 environment
```bash
# Activate qiime2 conda environment
source activate q2test

# Install plugin from github
pip install https://github.com/twbattaglia/q2-summarizer/zipball/master.zip

# Check plugin was installed
qiime info

# Check plugin is working correctly
qiime summarizer --help
```

## Usage
Download and convert example OTU-table
```
# Download example (.biom) file
curl -sO https://dl.dropboxusercontent.com/u/2868868/data/qiime2/feature-table.biom

# Convert biom table to QIIME2 Artifact
qiime tools import --type "FeatureTable[Frequency]" --input-path feature-table.biom --output-path feature-table-frequency.qza
```

Run summarizer on OTU-table to generate an interactive report
```bash
# Run summarizer on an otu-table
qiime summarizer otu-table --i-table feature-table-frequency.qza --o-visualization summarize_out

# View results
qiime tools view summarize_out.qzv
```

## History

TODO: Write history

## Credits
- Qiime2 https://github.com/qiime2/qiime2
- AdminLTE https://github.com/almasaeed2010/AdminLTE
