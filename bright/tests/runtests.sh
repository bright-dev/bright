#!/bin/bash
time for x in $(ls test_*.py); do nosetests $x; done
