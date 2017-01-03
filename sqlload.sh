#!/bin/bash

gunzip -c data/sqldump.txt.gz | mysql -u djangouser -p djtest
