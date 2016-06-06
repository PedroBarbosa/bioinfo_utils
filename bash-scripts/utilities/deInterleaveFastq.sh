#!/bin/bash
# Set up some defaults
GZIP_OUTPUT=0
PIGZ_COMPRESSION_THREADS=10

# If the third argument is the word "compress" then we'll compress the output using pigz
if [[ $3 == "compress" ]]; then
  GZIP_OUTPUT=1
fi

if [[ ${GZIP_OUTPUT} == 0 ]]; then
  paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > $1) | cut -f 5-8 | tr "\t" "\n" > $2
else
  paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $1) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $2
fi