#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# fetch the CLUES code and checkout the PALM branch
cd ..
git clone git@github.com:35ajstern/clues.git
cd clues/
git checkout aaron/palm-integration