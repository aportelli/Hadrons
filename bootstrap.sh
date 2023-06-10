#!/usr/bin/env bash

set -euo pipefail

sqlite_link="https://www.sqlite.org/2023/sqlite-amalgamation-3410200.zip"
sha256='01df06a84803c1ab4d62c64e995b151b2dbcf5dbc93bbc5eee213cb18225d987'

echo '-- generating module list...'
cd Hadrons
./make_module_list.sh
cd ..
echo '-- downloading SQLite...'
wget ${sqlite_link}
archive="$(basename ${sqlite_link})"
folder="${archive%.*}"
echo "${sha256} ${archive}" | sha256sum --check
unzip "${archive}"
mkdir -p Hadrons/sqlite/
mv "${folder}"/sqlite3* Hadrons/sqlite/
rm -rf "${folder}" "${archive}"
echo '-- generating configure script...'
mkdir -p .buildutils/m4
autoreconf -fvi
