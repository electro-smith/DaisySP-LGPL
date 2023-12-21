#!/bin/bash

# Check for the required project name
if [[ $# -lt 1 || ${1:0:1} == "-" ]]; then
  echo "Usage: $0 PROJECT_NAME [-l libdaisy_location] [-d daisysp_location]" >&2
  exit 1
fi

# Set the project name value
PROJECT="$1"

shift

# Parse command line options
while [[ $# -gt 0 ]]; do
  case $1 in
    -l)
      LIBDAISY_DIR="$2"
      shift 2
      ;;
    -d)
      DAISYSP_DIR="$2"
      shift 2
      ;;
    *)
      echo "Usage: $0 PROJECT_NAME [-l libdaisy_location] [-d daisysp_location]" >&2
      exit;
      ;;
  esac
done

# Default variable value
: "${LIBDAISY_DIR:="../../../libDaisy"}"
: "${DAISYSP_DIR:="../../../DaisySP"}"

# Make the folders if they don't exist
mkdir -p distribution/
mkdir -p distribution/licenses
mkdir -p distribution/resource

# copy .a files
cp $DAISYSP_DIR/build/libdaisysp.a distribution/resource || { echo 'Invalid libDaisy location. Try -l flag' ; exit 1; }
cp $DAISYSP_DIR/DaisySP-LGPL/build/libdaisysp-lgpl.a distribution/resource || { echo 'Invalid DaisySP location. Try -d flag' ; exit 1; }
cp $LIBDAISY_DIR/build/libdaisy.a distribution/resource || { echo 'Invalid DaisySP location. Try -d flag' ; exit 1; }

# copy bootloader files
cp $LIBDAISY_DIR/core/dsy_bootloader* distribution/resource || { echo 'Invalid libDaisy location. Try -l flag' ; exit 1; }
cp $LIBDAISY_DIR/core/STM32H750IB_flash.lds distribution/resource || { echo 'Invalid libDaisy location. Try -l flag' ; exit 1; }
cp $LIBDAISY_DIR/core/STM32H750IB_sram.lds distribution/resource || { echo 'Invalid libDaisy location. Try -l flag' ; exit 1; }
cp $LIBDAISY_DIR/core/STM32H750IB_qspi.lds distribution/resource || { echo 'Invalid libDaisy location. Try -l flag' ; exit 1; }

# copy and prepend vars to local makefile
cp $DAISYSP_DIR/DaisySP-LGPL/distribution/Makefile_LGPL distribution/Makefile || { echo 'Invalid libDaisy location. Try -l flag' ; exit 1; }

prepend="# Project Name\nTARGET = $PROJECT\n"
sed -i -e "1i $prepend" distribution/Makefile

# copy .o and .map files
cp build/*.o distribution/resource
cp build/*.map distribution/resource

# Copy LICENSE files
cp $DAISYSP_DIR/DaisySP-LGPL/LICENSE distribution/licenses/DaisySP-LGPL
cp $DAISYSP_DIR/LICENSE distribution/licenses/DaisySP
cp $LIBDAISY_DIR/LICENSE distribution/licenses/libDaisy

# make a little readme
readme="To rebuild with an updated or modified version of any of the static libraries:\n
1. Install the Daisy-Toolchain.\n
2. Build the library you would like to replace, for example by running make inside of DaisySP-LGPL.\n
3. Copy the generated .a file to the resources folder. In this case it would be libdaisy-lgpl.a \n
4. Run the make command in this folder. Your .bin file will be in the build folder. "

echo -e $readme > distribution/Readme.md