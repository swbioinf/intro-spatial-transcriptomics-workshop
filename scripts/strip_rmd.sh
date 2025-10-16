#!/bin/bash -ue

input_file=$1
out="notebooks/${input_file}"

mkdir -p notebooks

sed -n '
	# Print code chunks, do not incldue comments again later
	/^```/,/```$/ {
		p
		/```$/a

		b
	}

	# Print headers, excluding comments above
	/^#/ {
		p
		a

	}
' "$1" > ${out}

echo "${input_file} → ${out}"
