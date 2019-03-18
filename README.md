# reticulatus
A long snake for a long assemblies

## How to drive this thing

#### Setup the environment

```
conda env create --name reticulatus --file environments/base.yaml
conda activate reticulatus
```

#### Engage the pipeline

Run the pipeline with `snakemake`, you **must** specify `--use-conda` to ensure that
any tools that require a special jail (*e.g.* for `python2`) are run far, far away
from everything else.
Set `j` to the highest number of processes that you can fill with snakes before
your computer falls over.

```
snakemake -j <available_threads> --reason --use-conda
```
