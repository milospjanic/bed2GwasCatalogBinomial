# bed2GwasCatalogBinomial

The P-values will be computed using binomial cumulative distribution function b(x;n,p) in R (dbinom function). We set the parameter n equal to the total number of GWAS SNPs in a particular GWAS phenotype. Parameter x was set to the number of GWAS SNPs for a given GWAS phenotype that overlap input regions and parameter p was set to the fraction of the uniquely mappable human hg19 genome (calculated with subscript) that is localized in the input regions. Calculated binomial p-value equals the probability of having x or more of the n test genomic regions in the input domains given that the probability of that occurring for a single genomic location is p. 

This script will connect to GWAS Catalog and download the entire data set, create bed file, and parse and uniq according to the N-1 input arguments. Last argument provided to the bash script should be a bed file that will be used to intersect parsed bed files from the GWAS Catalog. Number of overlaps is reported and initial number of entries in parsed files.  Finally, it will create an R script that will be executed to calculate binomial p-values for each overlap. Intermediary files will be removed, except: GwasCatalog.bed (entire catalog in a bed file), \*gwascatalog.bed (parsed original files from GWAS Catalog), \*gwascatalog.bed.cut.sort.uniq.chrXY (parsed and uniqed files from GWAS Catalog), \*gwascatalog.bed.cut.sort.uniq.overlap (overlap with parsed files), GwasCatalog2Bed.sh (sciprt to download GWAS Catalog and convert to bed). Note that names of phenotypes in GWAS Catalog start with capital letter but then next word is with small letter. **That is why I enabled case insensitive search in this script.**

#Usage
<pre>
chmod 775 ./bed2GwasCatalogBinomial.sh
./bed2GwasCatalogBinomial "Coronary_artery" "Coronary_heart" "Bipolar_disorder" "Feminism" L2_TCCTGAGC_L002_peaks.bed 
</pre>

#Dependencies 
Rscript, bedtools (needs to be in $PATH)

#Output
<pre>
./bed2GwasCatalogBinomial "Coronary_artery" "Coronary_heart" "Bipolar disorder" "Feminism" L2_TCCTGAGC_L002_peaks.bed 

--2016-04-20 14:41:34--  http://www.genome.gov/admin/gwascatalog.txt
Resolving www.genome.gov (www.genome.gov)... 156.40.242.24
Connecting to www.genome.gov (www.genome.gov)|156.40.242.24|:80... connected.
HTTP request sent, awaiting response... 301 Moved Permanently
Location: https://www.genome.gov/admin/gwascatalog.txt [following]
--2016-04-20 14:41:34--  https://www.genome.gov/admin/gwascatalog.txt
Connecting to www.genome.gov (www.genome.gov)|156.40.242.24|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 10407265 (9.9M) [text/plain]
Saving to: `gwascatalog.txt'

100%[=======================================================================================================================================================>] 10,407,265  1.49M/s   in 7.0s    

2016-04-20 14:41:42 (1.42 MB/s) - `gwascatalog.txt' saved [10407265/10407265]

Received: Coronary_artery
Done: Coronary_artery
Received: Coronary_heart
Done: Coronary_heart
Received: Bipolar_disorder
Done: Bipolar_disorder
Received: Feminism
Done: Feminism
Gwas Catalog number of SNP-phenotype associations:
18899 GwasCatalog.bed
Gwas Catalog number of SNP-phenotype associations per category:
Phenotype: Coronary_artery
110 Coronary_artery.gwascatalog.bed
Phenotype: Coronary_heart
143 Coronary_heart.gwascatalog.bed
Phenotype: Bipolar_disorder
444 Bipolar_disorder.gwascatalog.bed
Phenotype: Feminism
12 Feminism.gwascatalog.bed
Gwas Catalog number of SNP-phenotype associations per category AFTER REMOVING DUPLICATES:
Phenotype: Coronary_artery
89 Coronary_artery.gwascatalog.bed.cut.sort.uniq
Phenotype: Coronary_heart
130 Coronary_heart.gwascatalog.bed.cut.sort.uniq
Phenotype: Bipolar_disorder
411 Bipolar_disorder.gwascatalog.bed.cut.sort.uniq
Phenotype: Feminism
12 Feminism.gwascatalog.bed.cut.sort.uniq
Converting Phenotype: Coronary_artery
Converting Phenotype: Coronary_heart
Converting Phenotype: Bipolar_disorder
Converting Phenotype: Feminism
L2_TCCTGAGC_L002_peaks.bed
Input phenotypes:
Coronary_artery Coronary_heart Bipolar_disorder Feminism
Overlapping Phenotype SNPs with input bed: Coronary_artery
Overlapping Phenotype SNPs with input bed: Coronary_heart
Overlapping Phenotype SNPs with input bed: Bipolar_disorder
Overlapping Phenotype SNPs with input bed: Feminism
Number of Overlapping Phenotype SNPs with input bed: Coronary_artery
4 Coronary_artery.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Coronary_heart
7 Coronary_heart.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Bipolar_disorder
8 Bipolar disorder.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Feminism
0 Feminism.gwascatalog.bed.cut.sort.uniq.overlap
Number of Overlapping Phenotype SNPs with input bed: Coronary_artery
4
Number of Overlapping Phenotype SNPs with input bed: Coronary_heart
7
Number of Overlapping Phenotype SNPs with input bed: Bipolar_disorder
8
Number of Overlapping Phenotype SNPs with input bed: Feminism
0
--2016-04-20 14:41:44--  https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
Resolving genome.ucsc.edu (genome.ucsc.edu)... 128.114.119.134, 128.114.119.132, 128.114.119.131, ...
Connecting to genome.ucsc.edu (genome.ucsc.edu)|128.114.119.134|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1971 (1.9K) [text/plain]
Saving to: `hg19.chrom.sizes.18'

100%[=======================================================================================================================================================>] 1,971       --.-K/s   in 0s      

2016-04-20 14:41:44 (55.3 MB/s) - `hg19.chrom.sizes.18' saved [1971/1971]

Human Genome size version hg19: 3137161264
Coverage of BED file 91380849
Fraction of hg19 0.0291285
[1] "Coronary_artery"
[1] 0.1810647
[1] "Coronary_heart"
[1] 0.06673914
[1] "Bipolar_disorder"
[1] 0.04603295
[1] "Feminism"
[1] 0.7013601
</pre>
