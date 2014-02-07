#!/bin/bash
python ComputeRotationAngle.py "$1" "$2" >& RotationAngle`if [ -n "$2" ]; then echo "_$2"; fi`.log
python MakeList.py "$1" "$2"
