#!/bin/bash

mysqldump -u djangouser -p djtest | gzip -c > data/sqldump.txt.gz
