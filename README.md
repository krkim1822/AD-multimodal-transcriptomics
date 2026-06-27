# Multimodal_Transcriptomics_AD

### Multimodal Analysis of Transcriptional Regulation in the Dorsolateral Prefrontal Cortex Enhances Alzheimer’s Disease Risk Gene Discovery

![image](https://github.com/krkim1822/AD-multimodal-transcriptomics/blob/main/Multimodal%20analyses%20framework.png)

## Quantify multimodal transcriptomic traits from raw RNA-seq data
### Environment Setup
* Clone the [Pantry](https://github.com/PejLab/Pantry) repository and install dependencies using:
  ```
  git clone https://github.com/PejLab/Pantry.git

  cd Pantry/phenotyping
  conda env create -n pantry --file environment.yml
  conda activate pantry
  ```
  Detailed instructions for installation, dependencies, and required input file formats can be found in the [Pantry](https://github.com/PejLab/Pantry). 

### Configuration
The `config.yml` file is used to specify general parameters, reference genomes, input files, and the specific transcriptomic modalities to generate. The default configuration file provided in [Pantry/phenotyping/](https://github.com/PejLab/Pantry/tree/main/phenotyping) must be customized for your own dataset. 
* Example config.yml
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
This pipeline is executed via the Snakemake workflow management system. Activate the environment, and run Snakemake inside the phenotyping directory:
* Example
  ```
  cd Pantry/phenotyping
  ### Ensure your configuration file is saved as 'config.yml' in this directory. 
  ### Run Snakemake using the specified parallel cores 
  snakemake --jobs 32 
  ```
