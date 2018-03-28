#!/bin/bash

# Run it as $bash cluster.sh 2.0 ../native.pdb ../decoys.silent

total_decoys=`wc -l ../score.sc | awk '{print $1}'`
clustering_pool_num=$((total_decoys/100))

printScoreFile_byHeader.pl ../score.sc I_sc reweighted_sc score_lowres_opt description | sort -nk 3 | head -n $clustering_pool_num | awk '{print $NF}' >pdb_list

if [ -f list ];then rm list;fi


PATH_TO_EXE="/vol/ek/Home/alisa/rosetta/Rosetta/main/source/bin"
PATH_TO_DB="/vol/ek/Home/alisa/rosetta/Rosetta/main/database/"

radius=$1
native=$2
decoys=$3

len=`grep CA $native | wc -l`
plen=`awk '$5=="B"' $native | grep CA | wc -l`
actualR=`date | awk '{print sqrt('$plen'/'$len')*'$radius'}'`
echo actual radius is $actualR

if [ ! -e cluster.silent ]; then
		$PATH_TO_EXE/cluster.linuxgccrelease -in:file:silent $decoys -in:file:silent_struct_type binary -database $PATH_TO_DB -cluster:radius $actualR -in:file:fullatom -tags `cat pdb_list` -silent_read_through_errors > clog
    echo Done clustering.
fi

echo Now printing results.

#x=`wc -l pdb_list | awk '{print $1}'`
#tail -$((x+5)) clog | head -${x} | awk '{print $4,$5,$6}' > cluster_list

x=`wc -l pdb_list | awk '{print $1}'`
tail -$((x+5)) clog | head -${x} | grep _0001 | awk  '{for(i=1;i<=NF;i++){if ($i ~ /_0001/){print $i,$(i+1),$(i+2)}}}' | sed '/^$/d' > cluster_list


if [ -e pdb_list_sc ];then rm pdb_list_sc;fi
for i in `awk '{print $1}' cluster_list`;do sed 1d ../score.sc | head -1 > tmp1; grep $i ../score.sc >> tmp1; printScoreFile_byHeader.pl tmp1 I_sc reweighted_sc score_lowres_opt | tail -1 | awk '{print $2,$3,$4}' >>pdb_list_sc; done
paste cluster_list pdb_list_sc >cluster_list_sc

echo "Decoy_ID Cluster_no Member_ID I_sc reweighted_sc" >cluster_list_I_sc_sorted
echo "Decoy_ID Cluster_no Member_ID I_sc reweighted_sc" >cluster_list_reweighted_sc_sorted

sort -nk 4 cluster_list_sc | sort -u -k2,2 | sort -nk 4 | head -20 >>cluster_list_I_sc_sorted
sort -nk 5 cluster_list_sc | sort -u -k2,2 | sort -nk 5 | head -20 >>cluster_list_reweighted_sc_sorted

exit 0
