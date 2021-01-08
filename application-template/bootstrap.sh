#!/usr/bin/env bash
 set -e
 
 echo '-- generating configure script...'
 mkdir -p .buildutils/m4
 autoreconf -fvi
 