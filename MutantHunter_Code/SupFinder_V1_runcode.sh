cd /Volumes/MITCH/Martin_fastq/

chmod 777 Supfinder_V1.sh

./Supfinder_V1.sh -m 14_S14 -n 13_S13 -g genome -f genome -r paired -s 50




cd /Volumes/MITCH/Martin_fastq/

chmod 777 Supfinder_V1_noloop.sh

./Supfinder_V1_noloop.sh -m 14_S14 -n 13_S13 -g genome -f genome -r paired -s 100



cd /Volumes/MITCH/Martin_fastq/

chmod 777 Supfinder_V1_noloop.sh

./Supfinder_V1_partial_loop.sh -m 14_S14 -n 13_S13 -g genome -f genome -r paired -s 10
