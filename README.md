# SNP selection workflow for low-coverage whole genome sequencing data for GTseq applications

The workflow described here is based on previous work detailed in [Development of SNP Panels from Low-Coverage Whole Genome Sequencing (lcWGS) to Support Indigenous Fisheries for Three Salmonid Species in Northern Canada](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.14040).
It features scripts originally developped by Xavier Dallaire and Claire Mérot [`angsd_pipeline`](https://github.com/xav9536/angsd_pipeline/tree/master).
The python scripts for SNP selection were built by Eric Normandeau (https://github.com/enormandeau/gtseq_scripts).

## Pipeline overview

1. Call 2 sets of SNPs
    * A set of SNPs using various filtering criteria (minor allele frequency, number of supporting reads, ...)
    * A set of **background** SNPs using less stringent filtering criteria. This one will be used later in the pipeline to remove candidate loci that are located in regions bearing multiple potential SNPs, as these may introduce noise.
2. Determine population structure and differentiation from PCA, admixture and Fst
3. Select informative SNPs 
    * Extract most differentiated SNPs across populations from pairwise allele frequency differences
    * Score SNPs based on various criteria (GC contents, number of neighbor background SNPs within 100 bp, ...) and filter sites
    * Remove pre-selected sites that map to more than one region of the genome
4. Validate selected SNP set with PCA and admixture analysis
    * Confirm that population structure inferred from selected SNPs corresponds to population structured inferred from the original SNP set
    * Confirm that population structure inferred from **random subsets** of the selected SNP set corresponds to population structured inferred from the original SNP set


## Prerequisites

### Files
* A list of chromosomes in `02_infos/chrs.txt`, used for parallelizing the SNP calling step across chromosomes
* Bed files of regions to exclude from calling
	* Sex-linked regions
	* Repeated regions: this list produced by first running `RepeatModeler` to detect repeated elements in the reference genome, then running `RepeatMasker` with the custom repeat library and converting the .gff output to bed (see scripts `RepeatModeler.sh` and `RepeatMasker.sh` in `01_scripts/utils`)
	* List of non-ATCG sites in the reference genome, which can be made with `seqkit locate -P -i -r -G -p "[^A^T^C^G]+" --bed 03_genome/genome.fasta | cut -f1-3 > 02_infos/regions_to_exclude/ambiguous_bases.bed`
* A reference genome (.fasta)
* A tab-seperated file listing sample ID and their assigned population (`02_infos/ID_pop.txt`)
* A list of populations included in the study (can be made by `less 02_infos/ID_pop.txt | cut -f2 | sort | uniq > 02_infos/pops.txt`)
* A list of bam file paths to mapped short reads for each sample (`02_infos/bam.filelist`). See the [wgs_sample_preparation workflow](https://github.com/enormandeau/wgs_sample_preparation) for detailed read mapping steps.

### Software
* `angsd 0.937`
* `pcangsd 1.10`
* `samtools 1.15`
* `R` >4.2
* `bedtools 2.31.1`
* `ngsParalog 1.3.4`
* `ngsadmix` v33 (https://raw.githubusercontent.com/ANGSD/angsd/master/misc/ngsadmix32.cpp)
* `python 3.7`
* `ncbiblast 2.6.0`
* `bedtools 2.31.1`

## Detailed walkthrough
1. Prepare required files for SNP calling: `01_prepare_regions.sh`
2. Call 2 sets of SNPs (both series of scripts can be ran at the same time or not - only the first SNP set is required for steps that preceed SNP selection)
    * Call SNP sites that respect filtering criteria across all chromosomes (`02_call_SNPs_all.sh`), filter for canonical and deviant sites (`03_call_SNPs_canonical.sh`), then concatenate sites from all chromosomes into a single file (`04_concat_SNPs_canonical.sh`).
    * Call background SNP sites using less stringent criteria, using scripts `02.1_call_SNPs_background.sh`, `03.1_call_SNPs_canonical_background.sh` and `04.1_concat_SNPs_canonical_background.sh`.
3. Perform PCA (`05_pca.sh`) and admixture analysis (`06_ngs_admix.sh`).
4. Compute minor allele frequencies (`07_maf_by_pop.sh`) and site allele frequency likelihoods (`07_saf_by_pop.sh`) in each population for calculating Fst (`08_fst_by_group.sh`).
5. Preselect SNPs: `09_preselect_SNPs.sh`
    * After running this script, look at plots generated in the `10_SNP_selection` folder to determine appropriate thresholds for min and max GC contents, max MAF, max number of neighbor background SNPs and max complexity, and adjust variables accordingly in `10_filter_map_SNPs.sh`, `11_select_SNPs.sh` and `12_validate_SNPs.sh`.
6. Filer preselected sites and map their flanking sequence against the reference to keep sites that map to only one region with at least 90% identity and over at least 160 bp: `10_filter_map_SNPs.sh`
7. Extract the final SNPs (`11_select_SNPs.sh`) by keeping only sites seperated by at least 200kb (`MIN_DIST`)
    * The `TARGET_SUM` variable is a arbitrary number that will determine the number of selected SNPs. Play around with this values until you get close to the target number of SNP for the panel. 
8. Validate the selected SNP set by doing PCAs and admixture analysis on the final SNP set itself, as well as on 10 random subsets of 100 SNPs from the final SNP set (`12_validation.sh`).  
