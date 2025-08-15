#!/bin/awk -f

BEGIN { OFS="," }  # Set the output field separator to a comma for CSV
{
    if (NR % 4 == 1) {  # Header line (starts with @)
        id = substr($0, 2) # Remove the leading '@'
    } else if (NR % 4 == 2) { # Sequence line
        sequence = $0
    } else if (NR % 4 == 0) { # Quality score line
        print id, sequence
    }
}
