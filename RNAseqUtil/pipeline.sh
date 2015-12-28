#!/bin/bash

nohup ./pipeline.py < pipeline.txt > pipeline_log 2>&1 &
