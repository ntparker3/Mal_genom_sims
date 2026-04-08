This is a very barebones Nextflow workflow for running Dcifer on a lot of 
simulation outputs. Two example inputs are included in `data` for demonstration 
purposes. These were copied from 
[PGEcore](https://github.com/PlasmoGenEpi/PGEcore).

The workflow can be run in full mode (runs all inputs) or test mode (runs first 
input).

# Usage

To run locally:

```
nextflow run main.nf -profile local
```

To run locally with Docker (requires Docker installed and the image 
pulled/built):

```
nextflow run main.nf -profile local,docker
```

To run locally in test mode:

```
nextflow run main.nf -profile local --mode test
```

To run on Rockfish in test mode:

```
sbatch launch_workflow_test.sh
```

To run on Rockfish in full mode:

```
sbatch launch_workflow_full.sh
```

# Repo Structure

Brief orientation to the files and directories in this repo:

## Files

- main.nf: This is the Nextflow file that defines the workflow.
- nextflow.config: This file primarily defines various profiles. Current options 
  are local, rockfish, docker, and singularity. local and rockfish define 
  process resource configurations for running in different environments. docker 
  and singularity configure Nextflow to run in the corresponding container.
- Dockerfile: Used to build a Docker image, from which you can build a 
  Singularity container, which is what you can actually run from on Rockfish. 
  **Note:** For some reason I can't get this to build, and I currently have the 
  workflow pointing to the container I use for GenE8 relatedness analyses. That 
  should work, but I can try to troubleshoot this more. The main reason that 
  might be important is that Dcifer is fixed at 1.3.1 in that other container 
  because I was having a weird issue after some updates Inna made. If you need a 
  newer Dcifer, we need to address this.
- launch\_workflow\_{mode}.sh: These two shell scripts are used for running the 
  workflow on Rockfish.

## Directories

- data: Default directory for inputs, but this is overriden for running on 
  Rockfish.
- results: Default directory for outputs, but this is overriden for running on 
  Rockfish.
- bin: Directory for scripts and notebooks needed by the workflow.
- work: This will be created once you run it once. This is where Nextflow puts 
  intermediates.
- slurm: Created on Rockfish to store Slurm logs. Look at the one matching your 
  Slurm job ID (for the main workflow) to see how things went.
