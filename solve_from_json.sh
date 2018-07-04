./colour_order -o-5 <(python kep_h_with_multiword_weights/kep_h.py -f $1 \
	-c effective:size:inverse3way:backarc:weight -s 1:1:1:1 -e 3 -n 2 -i)
