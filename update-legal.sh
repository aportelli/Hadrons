#!/bin/bash

update_files() 
{
    while (( "$#" )); do

    echo $1

    cat > message  <<EOF
/*
 * $(basename $1), part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
EOF

    git log $1 | grep Author > gitauth
    grep 'Author: '  $1 > fileauth

    cat gitauth fileauth | awk '{if ($1 != "*"){printf(" * ")}; print $0}' | sort -u >> message

    rm gitauth fileauth

    cat >> message <<EOF
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
EOF

    cat message > tmp.fil

    NOTICE=`grep -n "END LEGAL" $1 | awk '{ print $1 }'  `

    if [ "X$NOTICE" != "X" ]
    then
        echo "found notice ending on line $NOTICE"
        awk 'BEGIN { P=0 } { if ( P ) print } /END LEGAL/{P=1} ' $1 >> tmp.fil
    else
        cat $1 >> tmp.fil
        
    fi


    cp tmp.fil $1

    shift

    done

    rm message tmp.fil
}

update_files $(find . -name '*.cpp')
update_files $(find . -name '*.hpp' | grep -v 'Modules.hpp')
