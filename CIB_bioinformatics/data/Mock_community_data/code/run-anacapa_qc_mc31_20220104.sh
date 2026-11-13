# EDIT THESE
BASEDIR="/data/anacapa" # change to folder you want shared into container
CONTAINER="/data/anacapa/anacapa-1.5.0.img" # change to full container .img path
DB="/data/anacapa/Anacapa-New-Master_08152019/Anacapa_db" # change to full path to Anacapa_db
DATA="/data/home/zgold/mock_community/Mockcommunity_20211231_v2" # change to input data folder (default 12S_test_data inside Anacapa_db)
OUT="/data/home/zgold/mock_community/Mockcommunity_20211231_output_v2" # change to output data folder

# OPTIONAL
FORWARD="/data/home/zgold/mock_community/forward_primers_fishcard_solo_p16s.txt"
REVERSE="/data/home/zgold/mock_community/reverse_primers_fishcard_solo_p16S.txt"
LENGTH="/data/home/zgold/mock_community/metabarcode_loci_min_merge_length_fishcard_solo_p16S.txt"
cd $BASEDIR

# If you need additional folders shared into the container, add additional -B arguments below

time singularity exec -B $BASEDIR $CONTAINER /bin/bash -c "$DB/anacapa_QC_dada2.sh -i $DATA -o $OUT -d $DB -f $FORWARD -r $REVERSE -e $LENGTH -a nextera -t MiSeq -l -q 30"