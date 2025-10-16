#!/bin/bash -ue

input_file=$1
out="notebooks/${input_file}"

mkdir -p notebooks

sed -n '
	# Print headers
	/^#/ {
		p
		a\

	}

	# Print code chunks
	/^```/,/```$/p

	/```$/ {
		a\

	}
' "$1" > ${out}

echo "${input_file} → ${out}"
