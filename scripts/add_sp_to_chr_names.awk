#!/bin/awk -f

BEGIN{
ARGC=2
if(ARGV[2] == "mouse"){prefix = ">Mouse_"}

else if(ARGV[2] == "human"){prefix = ">Human_"}

else exit 1
}

/>/{sub(/>/, prefix)}1


