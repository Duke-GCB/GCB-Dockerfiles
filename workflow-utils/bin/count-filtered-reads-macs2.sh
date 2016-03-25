#!/bin/bash

grep -m 1 "tags after filtering in treatment" $1 | awk '{print $NF}'