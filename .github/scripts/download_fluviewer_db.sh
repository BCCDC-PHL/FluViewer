#!/bin/bash

database_version="v0.1.8"

if [ -f .github/data/fluviewer_db-${database_version}/FluViewer_db.fa ]; then
    rm .github/data/fluviewer_db-${database_version}/FluViewer_db.fa
fi

mkdir -p .github/data/fluviewer_db-${database_version}

wget -O .github/data/fluviewer_db-${database_version}/FluViewer_db.fa.gz \
     https://github.com/BCCDC-PHL/FluViewer-db/raw/${database_version}/FluViewer_db.fa.gz

gunzip  .github/data/fluviewer_db-${database_version}/FluViewer_db.fa.gz
