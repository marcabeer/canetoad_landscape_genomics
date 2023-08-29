#!/bin/bash

ALLSNPS=ct_bayenv_05022022.txt
mkdir -p tempdir
split -a 10 -l 2 $ALLSNPS ./tempdir/snp_batch
ls tempdir > split_file_name.txt
