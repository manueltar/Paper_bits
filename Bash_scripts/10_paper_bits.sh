#!/bin/bash

analysis=$1
MASTER_ROUTE=$2

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

rm -rf $output_dir
mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files


#### enformer_graphs #############################


type=$(echo "enformer_graphs""_""$analysis")
outfile_enformer_graphs=$(echo "$Log_files""outfile_""$type"".log")
touch $outfile_enformer_graphs
echo -n "" > $outfile_enformer_graphs
name_enformer_graphs=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_enformer_graphs=$(echo "$Rscripts_path""34_MPRA_summary_2.R")

MPRA_results=$(echo "/group/soranzo/manuel.tardaguila/Enformer/Table_S5_MPRA_Results.tsv")


myjobid_enformer_graphs=$(sbatch --job-name=$name_enformer_graphs --output=$outfile_enformer_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_enformer_graphs --MPRA_results $MPRA_results --type $type --out $output_dir")
myjobid_seff_enformer_graphs=$(sbatch --dependency=afterany:$myjobid_enformer_graphs --open-mode=append --output=$outfile_enformer_graphs --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_enformer_graphs >> $outfile_enformer_graphs")
