# Multimodal_Transcriptomics_AD

### Multimodal Analysis of Transcriptional Regulation in the Dorsolateral Prefrontal Cortex Enhances Alzheimer’s Disease Risk Gene Discovery

![image](https://github.com/krkim1822/AD-multimodal-transcriptomics/blob/main/Multimodal%20analyses%20framework.png)

---

## Quantify multimodal transcriptomic traits from raw RNA-seq data
### Environment setup
* Clone the [Pantry](https://github.com/PejLab/Pantry) repository and install dependencies using:
  ```sh
  git clone https://github.com/PejLab/Pantry.git

  cd Pantry/phenotyping
  conda env create -n pantry --file environment.yml
  conda activate pantry
  ```
  Detailed instructions for installation, dependencies, and required input file formats can be found in the [Pantry](https://github.com/PejLab/Pantry). 

### Configuration
* The `config.yml` file is used to specify general parameters, reference genomes, input files, and the specific transcriptomic modalities to generate. The default configuration file provided in [Pantry/phenotyping/](https://github.com/PejLab/Pantry/tree/main/phenotyping) must be customized for your own dataset. 
  * Example 
  ```YAML
  ## Raw RNA seq data
  paired_end: True
  read_length: 151
  fastq_dir: input/fastq
  fastq_map: input/fastq_map.txt

  ## Reference files
  ref_genome: input/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
  ref_anno: input/ref/Homo_sapiens.GRCh38.113.gtf

  ## Samples to include
  samples_file: input/samples.txt

  ## Modality groups
  ### List the modality groups to extract, with list of expected output files for each. 
  modality_groups:
    expression:
      files:
      - expression.bed.gz
      - expression.bed.gz.tbi
      - isoforms.bed.gz
      - isoforms.bed.gz.tbi
      - isoforms.phenotype_groups.txt
  ...
  ```

### Execution
* This pipeline is executed via the Snakemake workflow management system. Activate the environment, and run Snakemake inside the phenotyping directory:
  * Example
  ``` sh
  cd Pantry/phenotyping
  ### Ensure your configuration file is saved as 'config.yml' in this directory. 
  ### Run Snakemake using the specified parallel cores 
  snakemake --jobs 32 
  ```

---

## xQTL mapping & TWAS weight computation using pheast
* The [Pheast](https://github.com/PejLab/Pantry/tree/main/pheast) pipeline, provided within the Pantry, runs cis-xQTLs mapping via tensorQTL and TWAS weight computation via FUSION. It utilizes the quantified multimodal RNA traits generated from the previous quantification steps. 

### Environment setup
* Install dependencies using:
  ```sh
  cd Pantry/phenotyping
  conda env create -n pantry --file environment.yml
  conda activate pantry
  ```

### Configuration 
* the `config.yml` file in the `pheast` directory is used to specify input files directories, modalities to use in downstream analyses, and the specific analyses to run.
  * Example
    ```YAML
    ## Points to the directory containing Pantry output BED files:
    Phenotype_dir: input/phenotypes
    ## Prefix to genotypes matrix file in plink1 format (bed/bim/fam):
    Geno_prefix: input/rosmap_wgs
    ## Samples to include:
    Samples_file: input/samples.txt

    ## Modalities to analyze (must be present in the phenotype_dir):
    modalities:
      cross_modality:
        grouped: true
      expression:
        grouped: false
      isoforms:
        grouped: true

    ## Analyses to run:
    ### 'qtl' runs xQTL mapping and 'twas' runs FUSION weight computation.
    analyses:
      qtl:
        files:
        - '{modality}.cis_qtl.txt.gz'
        - '{modality}.cis_independent_qtl.txt.gz'
        - '{modality}.trans_qtl.txt.gz `
      twas:
        files:
        - '{modality}.tar.bz2'
    ```
* **Cross-Modality Mapping** : If you are running cross-modality xQTL mapping, concatenate all RNA trait tables into one, along with a groups file that specifies all RNA traits per gene as a single group. Check `Pantry/phenotyping/scripts/combine_modalities.sh` for more details.  

### Custom covariates for analyses
* We implemented a customized covariate setting in our study (e.g., 5 trait PCs, 5 genotype PCs, and additional sample covariates including ID, batch, RIN, study type, education, sex, and age at death), whereas the default `pheast` pipeline typically utilizes 20 trait PCs and 5 genotype PCs.
* `pheast/steps/covariates.smk`: Adjusted input requirements and PC thresholds.
  ``` python
  rule covariates:
    """Compute genotype and expression PCs and combine."""
    input:
        geno = multiext(str(interm_dir / 'covar' / 'geno_pruned'), '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        extra = 'input/rosmap_covar.tsv', # Points to additional sample covariates
    output:
        interm_dir / 'covar' / '{modality}.covar.tsv',
    params:
        pruned_prefix = interm_dir / 'covar' / 'geno_pruned',
        n_geno_pcs = 5, # Specify the number of genotype PCs
        n_pheno_pcs = 5, # Specify the number of phenotype PCs
    resources:
        mem_mb = 16000,
    shell:
        """
        # {input.extra} passes the additional sample covariates to the R script
        Rscript scripts/covariates.R \
            {params.pruned_prefix} \
            {input.bed} \
            {params.n_geno_pcs} \
            {params.n_pheno_pcs} \
            {input.extra} \
            {output}
        """
  ```
* `covariates.R`: Combined the trait/genotype PCs with the sample covariates to generate the final covariate data. The custom script can be found in the `scripts/` directory. 

* In addition to the input files required by pheast, the sample_covar.tsv file was also placed in the input directory. This is a tab-separated file where IDs match the RNA trait and genotype files.
  * Example of `input/sample_covar.tsv`
    | ID | batch | rin | study | educ | msex | age_death |
    | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
    | SAMPLE_001 | Batch_A_01 | 7.8 | ROS | 16 | 1 | 85.23 |
    | SAMPLE_002 | Batch_A_01 | 8.2 | MAP | 22 | 0 | 91.45 |
    | SAMPLE_003 | Batch_A_02 | 6.9 | ROS | 12 | 1 | 78.89 |

  * Modify the column names (e.g., batch, study) inside the `covariates.R` script to match your specific covaraite variables. 

### Execution
```sh
cd Pantry/pheast
snakemake --jobs 32
```

---



