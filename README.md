# reticulatus
**A snakemake-based pipeline for assembling and polishing long nanopore reads**

Reticulatus was developed in part to manage the execution of [long-read mock community experiments](https://lomanlab.github.io/mockcommunity/) at the Loman Lab.
It turns out that it's quite good, so I've generalised it for any long-read nanopore experiments, so you too can enjoy highly-contiguous, blisteringly fast, cutting-edge assembly and polishing too.
Reticulatus was designed for [assembly of metagenomic data](https://academic.oup.com/gigascience/article/8/5/giz043/5486468), but we have tried it on the [odd isolate too](https://twitter.com/samstudio8/status/1169293404943081473).

Reticulatus **is not** an assembler or polisher, but a [well stacked set of bioinformatics blocks](https://twitter.com/sienkieee/status/1192876481942294530). Reticulatus tries to codify what we at the Loman Lab think is the current best-practice for nanopore bioinformatics into a (hopefully) easy-to-use pipeline, taking advantage of all the goodness of Snakemake while adding a few features; including:

* a text-based read config that allows automated simple read pre-processing (deduplication, subsampling, merging)
* a text-based run config that provides a trivial way to define assembly and polishing strategies
* automatic generation of assembly bandage-art
* very fast GPU-accelerated polishing (racon, medaka)
* automated reporting of coverage and identity for contigs, for a set of references

As an attempted embodiment of best practice, Reticulatus is under development all of the time. Feel free to open an issue if it looks broken or send a pull request if it could work better.

Just so you know, the development of Reticulatus has:

* [helped make `racon` even faster](https://github.com/clara-genomics/racon-gpu/issues/3)
* [demonstrated GPU accelerated tools can work on ONT hardware](https://github.com/clara-genomics/racon-gpu/issues/2) and made [the containers to do so, freely available](https://github.com/SamStudio8/reticulatus-containers/)
* [led to a port of `minidot` that works with `minimap2`](https://github.com/SamStudio8/minidot)
* led to a more efficient implementation of BAM-based read subsampling in pomoxis
* pushed some minor fixes to Snakemake


## How to drive this thing

### (0) Clone the repository where you want the magic to happen

```
git clone https://github.com/SamStudio8/reticulatus.git; cd reticulatus;
```

### (1) Setup the environment

```
conda env create --name reticulatus --file environments/base.yaml
conda activate reticulatus
cp Snakefile-base Snakefile
```

You will almost certaintly want the `Snakefile-base` rule set for the time being. Run `Snakefile-zymo`to replicate our [mock community benchmarking pipeline](https://github.com/LomanLab/mockcommunity).

**Note** It is important that you ensure `snakemake-minimal` package is installed automatically using the environment specified above. Not only is this easier, but makes sure that the version installed is suitable for the overriden `shell.py` that ships with `reticulatus`.

### (2) Write your configuation

```
cp config.yaml.example config.yaml
```

Replace the YAML keys as appropriate. Keys are:

| Key | Type | Description | 
|-----|------|-------------|
| `dehumanizer_database_root` | Path, optional | empty directory in which to download the dehumanizer references (requires ~8.5GB), you can ignore this if you're not going to remove contigs assigned as human by `kraken2` |
| `kraken2_database_root` | Path | path to pre-built kraken2 database (*i.e.* the directory containing the `.k2d` files), or the path to a directory in which to `wget` a copy of our 30GB [microbial database](https://lomanlab.github.io/mockcommunity/mc_databases.html). **If the database already exists, you must** `touch k2db.ok` in this directory or **bad things** will happen |
| `slack_token` | str, optional | if you want to be bombarded with slack messages regarding the success and failure of your snakes, insert a suitable bot API token here |
| `slack_channel` | str, optional | if using a `slack_token`, enter the name of the channel to send messages, including the leading `#` |
| `cuda` | boolean | set to `False` if you do not want GPU-acceleration and `True` if you have the means to go very fast (*i.e.* you have a CUDA-compatible GPU) |
| `medaka_env` | URI | path to a singularity image (simg) or sandbox container to run medaka (GPU) |

### (3) Tell reticulatus about your reads

```
cp reads.cfg.example reads.cfg
```

For each sample you have, add a tab delimited line with the following fields: 

| Key | Type | Description | 
|-----|------|-------------|
| `sample_name` | str | a unique string that can be used to refer to this sample/readset later |
| `ont` | Path* | path to your long reads |
| `i0` | Path*, optional | path to your single-pair short reads for this sample, otherwise you can just set to `-` |
| `i1` | Path*, optional | path to your left paired-end short reads |
| `i2` | Path*, optional | path to your right paired-end short reads |
| `*` | - | an arbitrary delimiter that has no purpose |
||| feel free to add your own columns for metadata here, fill your boots, reticulatus doesn't care |

**\*** You can pre-process reads by modfying their file path as follows:

| Option | Syntax | Description | 
|--------|--------|-------------|
| Remove duplicates | myreads.**rmdup**.fq.gz | remove reads with a duplicate sequence header (to fix occasional duplicate reads arising from basecalling) |
| Subset reads | myreads.**subset-N**.fq.gz | select a random subsample of `N%` (with integer `N` between 1-99) |
| Merge reads | /path/to/merged/reads/:myreads.fq.gz,myotherreads.fq.gz,... | a root path for merged reads, followed by a colon and a comma delimited list of files to `cat` together, the filename will be chosen automatically and you should not be upset by this |

Pre-processing can be chained, for example: `myreads.rmdup.subset-25.fq.gz`, will remove sequence name duplicates and take 25% of the result. You may also use this syntax to pre-process files for merging. Reticulatus will work out what needs to be done to generate the new read files, and will only need to do so once; even when you run the pipeline again in the future.
The processed reads will be written to the same directory as the original reads. Once this has been done, you can delete the original reads yourself, if you'd like.

**Important** If you're using the GPU, you must ensure the directories that contain your reads are bound to the singularity container with `-B` in `--singularity-args`, use the same path for inside as outside to make things easier.


### (4) Tell reticulatus about your plans

```
cp manifest.cfg.example manifest.cfg
```

For each pipe you want to run, add a tab delimited line with the following fields:

| Key | Type | Description | 
|-----|------|-------------|
| `uuid` | str | a unique identifier, it can be anything, it will be used as a prefix for every file generated by this pipe, **do not insert the `.` character here if you want things to work** |
| `repolish` | str | if you wish to reuse an assembly for a different polishing scheme, enter the corresponding `uuid` name here, otherwise it **must** be set to `-` |
| `samplename` | str | the read set to assemble and polish, it must be a key from `reads.cfg` |
| `spell` | str | the "spell" to configure your assembly and polishing, corresponding to a named configuration in `spellbook.py` |
| `polishpipe` | str | a minilanguage that determines the polishing strategy. strategies are of the format `<program>-<readtype>-<iterations>` and are chained with the `.` character. *e.g.* `racon-ont-4.medaka-ont-1.pilon-ill-1` will perform four rounds of iterative `racon` long-read polishing, followed by one round of medaka long-read polishing and finally one round of `pilon` short-read polishing. Currently the following polishers are supported: racon, medaka, pilon and dehumanizer. No polishing can be acheived by setting to `-`. |
| `medakamodel` | str | the option to pass to `medaka_consensus -m`, this corresponds to the model to use for medaka long-read polishing, it will depend on your ONT basecaller version |
|       |               | feel free to add your own columns for metadata here, fill your boots, reticulatus doesn't care |
| `cpu` | int, optional | override the number of available CPU cores to this limit. this is optional, but if you use the field and don't want to override a sample, you **must** specify `-` |
| `gpu` | int, optional | override the number of available GPU interfaces to this limit. this is optional, but if you use the field and don't want to override a sample, you **must** specify `-` |


### (5) Engage the pipeline

Run the pipeline with `snakemake`, you **must** specify `--use-conda` to ensure that
any tools that require a special jail (*e.g.* for `python2`) are run far, far away
from everything else.
Set `j` to the highest number of processes that you can fill with snakes before
your computer falls over.

#### Simple

```
snakemake -j <available_threads> --reason
```

#### Advanced (GPU)

Additionally you **must** specify `--use-singularity` to use containers **and** provide suitable `--singularity-args` to use the GPU and bind directories. You must bind the directory into which you have cloned reticulatus, as well as any other directories that contain your reads. Set the `dir_inside` and `dir_outside` keys to the same path to ensure the file paths inside the container, match those on the outside of the container.

*e.g.* 
```
'--nv -B /data/sam-projects/reticulatus-testing/:/data/sam-projects/reticulatus-testing/ -B /path/to/reads/dir/:/path/to/reads/dir/ -B /path/to/more/reads/dir/:/path/to/more/reads/dir/'
```

You should also set `--resources gpu=N` where `N` is the number of GPU interfaces shown in `nvidia-smi`.
Don't forget to use the GPU, you must set the `cuda` key to True in `config.cfg`.

```
snakemake -j <available_threads> --reason --use-conda --use-singularity --singularity-args '--nv -B <dir_inside>:<dir_outside>' -k --restart-times 1 --resources gpu=N
```

Using the GPU will accelerate the following steps:

* `polish_racon`: you will need a racon binary compiled with `CUDA`, for your system
* `polish_medaka`: you will need to specify an appropriate singularity container, or install medaka with GPU support yourself


## Housekeeping

Unless otherwise stated by a suitable header, the files within this repository are made available under the MIT license. If you use this pipeline, an acknowledgement in your work would be nice... Don't forget to [cite Snakemake](https://snakemake.readthedocs.io/en/stable/project_info/citations.html).

## Support

If reticulatus has saved your computing bill, [maybe buy me a beer](https://www.buymeacoffee.com/samstudio8)?
