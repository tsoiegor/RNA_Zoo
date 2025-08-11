#!/bin/awk -f


/>/{
prefix=">"$2"_"$3" "
sub(/>/, prefix)
print $1"_"$2
}

$1 !~ />/{
print $0
}

