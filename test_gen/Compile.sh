#!/bin/bash
echo "Starting the script..."
date
cmake ./
make
./run
echo "Done!"
