#!/bin/bash

expr $(cat $1 | wc -l) / 4

