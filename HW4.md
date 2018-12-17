# Homework #4
Cynthia I. Rodriguez
***
# Summarize partitions of a genome assembly

>Go to the most current download genomes section and download the gzipped fasta 
file for all chromosomes.

### ANSWER:
``` 
wget “ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz”
```
### Calculate the following for all sequences ≤ 100kb and all sequences > 100kb:
1. Total number of nucleotides
2. Total number of Ns
3. Total number of sequences
### ANSWER:
```
#Use module jje/kent the faSize tool
module load jje/jjeutils jje/kent
```

> To get sequences ≤ 100kb:
```
faFilter -maxSize=100000 dmel-all-chromosome-r6.24.fasta \
> faSize dmel_fasta_less_100000
```
> To calculate number of nucleotides, Ns and number of sequences for ≤ 100kb:
```
faSize dmel_fasta_less_100000 

Results:
6178042 bases (662593 N's 5515449 real 5515449 upper 0 lower) in 1863 sequences in 1 files
Total size: mean 3316.2 sd 7116.2 min 544 (211000022279089) max 88768 (Unmapped_Scaffold_8_D1580_D1567) median 1567
N count: mean 355.7 sd 1700.6
U count: mean 2960.5 sd 6351.5
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real
```
Total number of nucleotides: 6178042, total number of N's: 662593, and total number of sequences: 1863.

> To get sequences >100kb:
faFilter -minSize=100000
```
dmel-all-chromosome-r6.24.fasta \
> faSize dmel_fasta_more_100000
```
> To calculate number of nucleotides, Ns and number of sequences for >100kb:
```
faSize dmel_fasta_more_100000

Results:
137547960 bases (490385 N's 137057575 real 137057575 upper 0 lower) in 7 sequences in 1 files
Total size: mean 19649708.6 sd 12099037.5 min 1348131 (4) max 32079331 (3R) median 23542271
N count: mean 70055.0 sd 92459.2
U count: mean 19579653.6 sd 12138278.9
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

```
Total number of nucleotides: 137547960, total number of N's: 490385, and total number of sequences: 7.

### Plots of the following for the whole genome, for all sequences ≤ 100kb, and all sequences > 100kb:
1. Sequence length distribution
2. Sequence GC% distribution
3. Cumulative genome size sorted from largest to smallest sequences
  ### 1. ANSWER FOR SEQUENCE LENGTH DISTRIBUTION:
 
 ##### For whole genome:
 ```
# Load the modules necessary for bioawk to work:
module load jje/jjeutils
module load jje/kent
bioawk -c fastx '{ print $name, length($seq) }' dmel-all-chromosome-r6.24.fasta > dmel_all_seqlenght.txt
```
Add the headers "Name" and "Length" to the txt file "dmel_all_seqlenght.txt" in excel then in RStudio do:
```
# To import the data into RStudio:
length_all$Cut <-cut(x=length_all$Length, breaks=10000)

#To load module:
library(ggplot2)

length_all_plot <-ggplot(data=length_all)

length_all_plot + geom_bar(mapping = aes(x=Cut)) + labs(title="Sequence length distribution for all sequences", x="Length", y="Number of Contigs") + theme_bw()+ theme(axis.text.x = element_text (angle = 90, hjust = 1))
```
IMAGE RESULT:
![Rplot_allseqlength.png](https://github.com/cirodri16/Homework-4/blob/master/Rplot_allseqlength.png?raw=true)

##### For all sequences ≤ 100kb:
```
bioawk -c fastx '{ print $name, length($seq) }' dmel_fasta_less_100000.fasta > dmel_less_seqlenght.txt
```
Add the headers "Name" and "Length" to the txt file "dmel_less_seqlenght.txt" in excel then in RStudio do:
```
length_less <- read.table("dmel_less_seqlenght.txt", header = TRUE)

length_less$Cut <-cut(x=length_less$Length, breaks=100)

library(ggplot2)

length_less_plot <-ggplot(data=length_less)

length_less_plot + geom_bar(mapping = aes(x=Cut)) + labs(title="Sequence length distribution for sequences less or equal 100kb", x="Length", y="Number of Contigs") + theme_bw()+ theme(axis.text.x = element_text (angle = 90, hjust = 1))
```
IMAGE RESULT:
![Rplot_less_seqlength.png](https://github.com/cirodri16/Homework-4/blob/master/Rplot_less_seqlength.png?raw=true)
##### For all sequences >100kb:
```
bioawk -c fastx '{ print $name, length($seq) }' dmel_fasta_more_100000.fasta > dmel_more_seqlenght.txt
```
Add the headers "Name" and "Length" to the txt file "dmel_more_seqlenght.txt" in excel then in RStudio do:
```
length_more <- read.table("dmel_more_seqlenght.txt", header = TRUE)

length_more$Cut <-cut(x=length_more$Length, breaks=100)

library(ggplot2)

length_more_plot <-ggplot(data=length_more)

length_more_plot + geom_bar(mapping = aes(x=Cut)) + labs(title="Sequence length distribution for sequences greater than 100kb", x="Length", y="Number of Contigs") + theme_bw()+ theme(axis.text.x = element_text (angle = 90, hjust = 1))
```
IMAGE RESULT:
![Rplot_more_seqlength.png](https://github.com/cirodri16/Homework-4/blob/master/Rplot_more_seqlength.png?raw=true)
### 2. ANSWER FOR SEQUENCE GC% DISTRIBUTION:
##### For whole genome:
```
bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.24.fasta
> dmel_all_GC.txt
```
Add headers to txt file in excel "Names" and "GC_Percent" then in RStudio:
```
GC_all <-read.table("dmel_all_GC.txt", header=TRUE)

GC_all$GC_Percentcut <-cut(x=GC_all$GC_Percent, breaks = 10)

library(ggplot2)

GC_all_plot <-ggplot(data=GC_all)

GC_all_plot + geom_bar(mapping = aes(x=GC_all$GC_Percentcut)) + labs(title="Sequence GC distribution for all sequences", x="GC Percentage", y="Number of Contigs") + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```
IMAGE RESULT:
![Rplot_allGC.png](https://github.com/cirodri16/Homework-4/blob/master/Rplot_all_GC.png?raw=true)

##### For all sequences ≤ 100kb:
```
bioawk -c fastx '{ print $name, gc($seq) }' dmel_fasta_less_100000.fasta
> dmel_less_GC.txt
```

Add headers to txt file in excel "Names" and "GC_Percent" then in RStudio:
```
GC_less <-read.table("dmel_less_GC.txt", header=TRUE)

GC_less$GC_Percentcut <-cut(x=GC_less$GC_Percent, breaks = 10)

library(ggplot2)

GC_less_plot <-ggplot(data=GC_less)

GC_less_plot + geom_bar(mapping = aes(x=GC_less$GC_Percentcut)) + labs(title="Sequence GC distribution for sequences equal or less than 100KB", x="GC Percentage", y="Number of Contigs") + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```
IMAGE RESULT:

![Rplot_lessGC.png](https://github.com/cirodri16/Homework-4/blob/master/Rplot_less_GC.png?raw=true)
##### For all sequences > 100kb:
```
bioawk -c fastx '{ print $name, gc($seq) }' dmel_fasta_more_100000.fasta
> dmel_more_GC.txt
```
Add headers to txt file in excel "Names" and "GC_Percent" then in RStudio:
```
GC_more <-read.table("dmel_more_GC.txt", header=TRUE)

GC_more$GC_Percentcut <-cut(x=GC_more$GC_Percent, breaks = 10)

library(ggplot2)

GC_more_plot <-ggplot(data=GC_more)

GC_more_plot + geom_bar(mapping = aes(x=GC_more$GC_Percentcut)) + labs(title="Sequence GC distribution for sequences greater than 100KB", x="GC Percentage", y="Number of Contigs") + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
```
IMAGE RESULT:
![Rplot_moreGC.png](https://github.com/cirodri16/Homework-4/blob/master/Rplot_moreGC.png?raw=true)



### 3. ANSWER FOR CUMULATIVE GENOME SIZE SORTED FROM LARGEST TO SMALLEST:
##### For whole genome:
```
bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nkblength_Ctg\t0" } { print "kblength_Ctg\t" $1 } ' \
>  seq_dmel_all.lengths
plotCDF2 seq_dmel_all.lengths seq_all.png
```
IMAGE RESULT:
https://github.com/cirodri16/Homework-4/blob/master/seq_all_cumulative.png
![seq_all_cumulative.png](https://github.com/cirodri16/Homework-4/blob/master/seq_all_cumulative.png?raw=true)
##### For all sequences ≤ 100kb:
```
bioawk -c fastx ' { print length($seq) } ' dmel_fasta_less_100000.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nkblength_Ctg\t0" } { print "kblength_Ctg\t" $1 } ' \
>  seq_dmel_less.lengths
plotCDF2 seq_dmel_less.lengths seq_less.png
```
IMAGE RESULT:
![seq_less_cumulative.png](https://github.com/cirodri16/Homework-4/blob/master/seq_less_cumulative.png?raw=true)

##### For all sequences >100kb :
```
bioawk -c fastx ' { print length($seq) } ' dmel_fasta_more_100000.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nkblength_Ctg\t0" } { print "kblength_Ctg\t" $1 } ' \
>  seq_dmel_more.lengths
plotCDF2 seq_dmel_more.lengths seq_more.png
```
IMAGE RESULT:
![seq_more_cumulative.png](https://github.com/cirodri16/Homework-4/blob/master/seq_more_cumulative.png?raw=true)

### Comments on "Summarize partitions of a genome assembly"

Perfection! As a matter of taste, you might consider using ```geom_histogram()```. It is more natural and easier to control for this purpose. But your plots are perfectly fine.

# Genome assembly
 > Assemble a genome from MinION reads:
1. Download the reads from here
2. Use minimap to overlap reads
3. Use miniasm to construct an assembly
### ANSWER:
```
module load jje/jjeutils
minimap=$(which minimap)
miniasm=$(which miniasm)
basedir=/pub/jje/ee282/$USER
projname=nanopore_assembly
basedir=$basedir/$projname
raw=$basedir/$projname/data/raw
processed=$basedir/$projname/data/processed
figures=$basedir/$projname/output/figures
reports=$basedir/$projname/output/reports

createProject $projname $basedir
ln -sf /bio/share/solarese/hw4/rawdata/iso1_onp_a2_1kb.fastq $raw/reads.fq

$minimap -t 32 -Sw5 -L100 -m0 $raw/reads.fq{,} \
| gzip -1 \
> $processed/onp.paf.gz

$miniasm -f $raw/reads.fq $processed/onp.paf.gz \
> $processed/reads.gfa

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| fold -w 60 \
> $processed/unitigs.fa
```
> Assembly assessment
1. Calculate the N50 of your assembly 
### ANSWER:
```
bioawk -c fastx \
 ' { li=length ($seq); l=li+l; print li; } END { print l; } ' \
unitigs.fa \
| sort -rn \
| gawk ' NR == 1 { l=$1; } 
NR >1 {
 li=$1; lc=li+lc; if(lc/l >= 0.5) {print li; exit; }
} ' \
| less -S
```
Result:
4494246

The above number is different to the Drosophila community reference's contig N50: 21,485,538 from NCBI.

2. Compare your assembly to the contig assembly (not the scaffold assembly!) from Drosophila melanogaster on FlyBase using a dotplot constructed with MUMmer (Hint: use faSplitByN as demonstrated in class)
### ANSWER:
```
# To split the scaffold assembly from ISO1 into contigs:

module load  jje/jjeutils perl
module load jje/kent/2014.02.19

faSplitByN dmel-all-chromosome-r6.24.fasta  contigs_dmel-all-chromosome-r6.24.fasta 10

# Check its N50 for splitting assembly into contigs:
bioawk -c fastx \ 
 ' { li=length ($seq); l=li+l; print li; } END { print l; } ' \
contigs_dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| gawk ' NR == 1 { l=$1; } 
NR >1 {
 li=$1; lc=li+lc; if(lc/l >= 0.5) {print li; exit; }
} ' \
| less -S
```
RESULT:
N50= 21485538
```
# Check its N50 for the original file:
bioawk -c fastx \ 
 ' { li=length ($seq); l=li+l; print li; } END { print l; } ' \
dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| gawk ' NR == 1 { l=$1; } 
NR >1 {
 li=$1; lc=li+lc; if(lc/l >= 0.5) {print li; exit; }
} ' \
| less -S
```
RESULT:
25286936


- Submitted the following script as a job in the HPC to compare my ONP assembly to the contig assembly that we just created:
```
#!/bin/bash
#$ -pe openmp 64
#$ -N test_HW4mumplot
#$ -q mic
#$ -m beas
#$ -M cirodri1@uci.edu
###Loading of binaries via module load or PATH reassignment
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot

###Query and Reference Assignment. State my prefix for output filenames
REF="contigs_dmel-all-chromosome-r6.24.fasta"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls u*.fa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)

###Used a value between 75-150 for -c. The value of 1000 is too strict.
nucmer -l 100 -c 75 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta \
  -R ${REF} -Q ${QRY} --postscript
```
Also submitted the job by changing the module load gnuplot/4.6.0 and the {QRY} --png and I obtained a png image.

To submit the job:
```
# Make file executable:
chmod u+x HW4codemum.bash

# Submit the job to the HPC cluster:
qsub ./HW4codemum.bash

# To see status of job:
qstat -u cirodri1

# To see the resulted image:
display flybase_unitigs.ps
OR
display flybase_unitigs.png
```
RESULT IMAGE:
ps image:
![mumer.PNG](https://www.dropbox.com/s/cg9fzctn08uficr/mumer.PNG?dl=0&raw=1)

PNG image:
![flybase2_unitigs.png](https://github.com/cirodri16/Homework-4/blob/master/flybase2_unitigs.png?raw=true)

3. Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot (Hint: use plotCDF2 as demonstrated in class and see this example)
### ANSWER:
```
module load rstudio/0.99.9.9
module load perl
module load jje/jjeutils/0.1a
module load jje/kent

# Oxford Nanopore assembly:
bioawk -c fastx ' { print length($seq) } ' /pub/jje/ee282/cirodri1/HW4/unitigs.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nUnitigs_Ctg\t0" } { print "Unitigs_Ctg\t" $1 } ' \
> unitigs.sorted.text

# Scaffold ISO1 from Flybase:
bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\ndmel_scaffold\t0" } { print "dmel_scaffold\t" $1 } ' \
> dmel_scaffold.sorted.text

# Contig assembly of ISO1 from Flybase:
bioawk -c fastx ' { print length($seq) } ' contigs_dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\ndmel_ctg\t0" } { print "dmel_ctg\t" $1 } ' \
> dmel_ctg.sorted.text

plotCDF2 *.sorted.text /dev/stdout \
| tee hw4_plotcdf2.png \
| display
```

RESULT IMAGE:
![CDF_assembly.png](https://github.com/cirodri16/Homework-4/blob/master/CDF_assembly.png?raw=true)

4. Calculate BUSCO scores of both assemblies and compare them- Submitted the following script as a job in the HPC:

### ANSWER:
```
#!/bin/bash
#$ -pe openmp 64
#$ -N test_HW4buscos
#$ -q abio,mic
#$ -m beas
#$ -M cirodri1@uci.edu

module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
QRY="unitigs.fa"
MYEXT=".fa" 

BUSCO.py -c 32 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}
```
To submit the job:
```
# Make file executable:
chmod u+x HW4codebuscos.bash

# Submit the job to the HPC cluster:
qsub ./HW4codebuscos.bash

# To see status of job:
qstat -u cirodri1
```
Results:
```
#To see results:
less short_summary_unitigs_diptera_odb9.txt

RESULTS:
# BUSCO version is: 2.0 
# The lineage dataset is: diptera_odb9 (Creation date: 2016-10-21, number of species: 25, number of BUSCOs: 2799)
# To reproduce this run: python /pub/jje/ee282/bin/busco/BUSCO.py -i unitigs.fa -o unitigs_diptera_odb9 -l /pub/jje/ee282/bin/busco/lineages/diptera_odb9/ -m genome -c 32 -sp fly
#
# Summarized benchmarking in BUSCO notation for file unitigs.fa
# BUSCO was run in mode: genome

        C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799

        13      Complete BUSCOs (C)
        13      Complete and single-copy BUSCOs (S)
        0       Complete and duplicated BUSCOs (D)
        32      Fragmented BUSCOs (F)
        2754    Missing BUSCOs (M)
        2799    Total BUSCO groups searched
```
The same script/code was submitted for the ISO1 from Flybase replacing the query file "QRY="unitigs.fa" to QRY="dmel-all-chromosome-r6.24"
Results:
```
# BUSCO version is: 2.0
# The lineage dataset is: diptera_odb9 (Creation date: 2016-10-21, number of species: 25, number of BUSCOs: 2799)
# To reproduce this run: python /pub/jje/ee282/bin/busco/BUSCO.py -i dmel-all-chromosome-r6.24.fasta -o dmel-all-chromosome-r6.24_diptera_odb9 -l /pub/jje/ee282/bin/busco/lineages/diptera_odb9/ -m genome -c 64 -sp fly
#
# Summarized benchmarking in BUSCO notation for file dmel-all-chromosome-r6.24.fasta
# BUSCO was run in mode: genome

        C:98.6%[S:98.1%,D:0.5%],F:0.8%,M:0.6%,n:2799

        2761    Complete BUSCOs (C)
        2747    Complete and single-copy BUSCOs (S)
        14      Complete and duplicated BUSCOs (D)
        21      Fragmented BUSCOs (F)
        17      Missing BUSCOs (M)
        2799    Total BUSCO groups searched
short_summary_dmel-all-chromosome-r6.24_diptera_odb9.txt (END)
```
The Flybase ISO1 assembly is more complete since it has a higher score for complete BUSCOs at 98.6% whereas the assembly with ONP data that we ran was at 0.5%.
