

grep -A 1 "snp_type_nonsynonymous" index.html > ./test.txt


grep -oP '(?<=align="center"><i>).*(?=</i>&nbsp;)' test.txt > genomiquegrep.txt



grep -A 1 "coding" index.html > ./insertdeletion.txt


grep -oP '(?<=align="center"><i>).*(?=</i>&nbsp;)' insertdeletion.txt > genomiquegrepinsertdelet.txt



grep -A 1 "snp_type_nonsense" index.html > ./stopcodon.txt

grep -oP '(?<=align="center"><i>).*(?=</i>&nbsp;)' stopcodon.txt > genomiquegrepstopcodon.txt

