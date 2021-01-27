#!/usr/bin/env bash

QP_VER=1.11.2
git clone https://github.com/arahlin/qpoint.git  --branch $QP_VER
cd qpoint
python3 setup.py install
