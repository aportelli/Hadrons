#!/usr/bin/env bash
set -e

echo '-- generating module list...'
cd Hadrons
./make_module_list.sh
cd ..
echo '-- generating configure script...'
mkdir -p .buildutils/m4
autoreconf -fvi
