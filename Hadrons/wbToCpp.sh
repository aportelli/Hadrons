#!/usr/bin/env bash

print_array () {
    local FILE=$1
    local li=$2
    local lf=$3
    awk -v li=${li} -v lf=${lf} 'NR == li {printf("    {%s,\n",$1)} NR > li && NR < lf {printf("    %s,\n",$1)} NR == lf {printf("    %s}",$1)}'  ${FILE}
}


if (( $# < 1 )); then
    echo "usage: `basename $0` <wavelet name1> [<wavelet name2> ...]"] 1>&2
    exit 1
fi

#tac ${LIST} | awk 'BEGIN{print "{"} {printf("%s,\n", $1)} END{print "}"}'

for NAME in $@; do
    TMP=`mktemp`
    curl -s http://wavelets.pybytes.com/wavelet/${NAME}/ | grep -Eo  -- '-?[0-9]+\.[0-9]+(e[+-][0-9]+)?$' | tac > ${TMP}
    size=$(( $(cat ${TMP} | wc -l)/4 ))
    echo "DwtFilter DwtFilters::${NAME} = {"
    echo "    // fwdL"
    print_array ${TMP} $(( 3*size + 1)) $(( 4*size )) 
    printf ',\n'
    echo "    // fwdH"
    print_array ${TMP} $(( 1*size + 1)) $(( 2*size ))
    printf ',\n'
    echo "    // bwdL"
    print_array ${TMP} $(( 2*size + 1)) $(( 3*size ))
    printf ',\n'
    echo "    // bwdH" 
    print_array ${TMP} 1 ${size}
    printf '\n};\n'
    echo ''
    rm -f ${TMP}
done
