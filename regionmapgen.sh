#!/bin/bash

#python3 regionmaptask.py forward
#mv ../regionmap.csv ../pg_regionmap_forward.csv
#mv ../regionmap.png ../pg_regionmap_forward.png
python3 regionmaptask.py reverse
mv ../regionmap.csv ../pg_regionmap_reverse.csv
mv ../regionmap.png ../pg_regionmap_reverse.png
