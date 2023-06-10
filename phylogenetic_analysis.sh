#phylogenetic analysis script

# README #
# The name of the BUSCO directory for each of the sample refers to the sample name
# seqtk, Gblocks, AMAS, MAFFT and IQTREE packages should be installed
# sample_names.txt input file should include the list of samples name (one name per one line)
# buscoids.txt input file should inclide the list of used orthologs (one ortholog ID per one line)
# name.lst input file should include the list of sequences to be inclided (one sequence header per), can be omitted if all sequences are used
# run script as "sh phylogenetic_analysis.sh WD", where "WD" is the directory the results should be stored to


#set directory
WD=$1
echo "The chosen directory is: "$WD""

# copy all needed files and label them
while IFS= read -r sample;
do
  echo $sample
  while IFS= read -r buscoid;
  do
    echo $buscoid
    for dir in ./$sample/run_arachnida_odb10/busco_sequences/*;
    do
      for file in $dir/$buscoid.faa;
      do
        cp $file $WD/${buscoid}_${sample}.faa
        sed -i "/^>/ s/$/_$sample/" $WD/${buscoid}_${sample}.faa
        #sed -i "1s/$/_$sample/" $WD/${buscoid}_${sample}.faa
      done
    done
  done < buscoinds.txt
done < sample_names.txt

# change working directory
cp buscoinds.txt $WD/buscoinds.txt
cp sample_names.txt $WD/sample_names.txt
cd $WD

#fix the issue of samples names duplication in header
while IFS= read -r sample;
do
  for file in *.faa;
  do
    sed -i "s/${sample}_${sample}/${sample}/g" $file
    echo $sample
  done  
done < sample_names.txt

#filter sequences by header according to name.lst file
for file in *.faa;
do
  seqtk subseq $file name.lst > "${file%%.*}".fasta
done

#fix sequences names, contig information to be deleted
for file in *.fasta;
do
  sed -i 's/NODE_[0-9]\+_length_[0-9]\+_cov_[0-9]\+\.[0-9]\+:[0-9]\+-[0-9]\+_//g' $file
  sed -i 's/[A-Z]\+[0-9]\+\.[0-9]\+:[0-9]\+-[0-9]\+_//g' $file
done

#catenate sequences of same orthologs
while IFS= read -r buscoid;
do
   cat ${buscoid}*.fasta >> ${buscoid}_cat.faa
done < buscoinds.txt

#local alignment for each of the ortholog
for file in *_cat.faa;
do
  mafft-linsi $file > "${file%%_cat*}"_aln.faa
done

#alignments' procession
for file in *_aln.faa;
do
  Gblocks $file -b5=h -e=.gb
done

#alignment's concatenation
AMAS.py concat -f fasta -d aa -i *.gb -u nexus --part-format nexus

#phylogenetic tree construction with ModelFinder tool, UFBootstrap 1000, SH-aLRT 1000, Nymphon striatum as outgroup, 2 threads used
iqtree -s concatenated.out -m MFP+MERGE -bb 1000 -alrt 1000 -rclusterf 10 -o Nymphon_striatum -nt 2 --prefix withpout_partition

#phylogenetic tree construction with Partition model, UFBootstrap 1000, SH-aLRT 1000, 2 threads used
iqtree -s concatenated.out -p partitions.txt -m MFP+MERGE -B 1000 -alrt 1000 -rclusterf 10 -nt 2 --prefix partition

#phylogenetic tree construction with C30 Mixture model, UFBootstrap 1000, SH-aLRT 1000, 2 threads used
iqtree -s concatenated.out -p partitions.txt -m LG+C30 -B 1000 -alrt 1000 -rclusterf 10 -nt 2 --prefix C30

#phylogenetic tree construction with C60 Mixture model, UFBootstrap 1000, SH-aLRT 1000, 2 threads used
iqtree -s concatenated.out -p partitions.txt -m LG+C60 -B 1000 -alrt 1000 -rclusterf 10 -nt 2 --prefix C60

#gCF/sCF analysis
iqtree -s concatenated.out -S partitions.txt --prefix loci -nt 2
iqtree -t partition.treefile --gcf loci.treefile -s concatenated.out --scf 100 --prefix concord

