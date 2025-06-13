import pandas as pd
import os
from glob import glob
from shutil import which
from snakemake.utils import min_version

### set minimum snakemake version ###
min_version("8.25.0")

meta = pd.read_csv("bin/genome_metadata.tsv", sep="\t")

rule all:
    input:
                # gtdbtk
                "analysis/gtdbtk/",
                # kraken2
                expand("analysis/kraken2/kraken2_{isolate}.out", isolate=meta["isolate"].tolist()),
                # kraken2_concat
                expand("analysis/kraken2/kraken2_concat_{isolate}.out", isolate=meta["isolate"].tolist()),
                # kraken2_summarize
                "analysis/kraken2/kraken2_summarize.txt",
                # kraken2_concat_summarize
                "analysis/kraken2/kraken2_concat_summarize.txt",

rule gtdbtk:
    input:      faa = "raw_data/fna/",
                batch = "bin/gtdbtk_batchfile.tsv",
    output:
                directory("analysis/gtdbtk/"),
    log:
                o = "logs/gtdbtk.o",
                e = "logs/gtdbtk.e",
    envmodules:
                "gtdbtk/2.4.0",
    resources:
                threads = 20,
                time = "48:00:00",
                nodes = 1,
                mem_gb = 128,
                partition = "priority-mem",
                qos = "pgec",
                name = "gtdbtk",
    shell:
                """
                gtdbtk classify_wf \
                  --genome_dir {input.faa} \
                  --skip_ani_screen \
                  --out_dir {output} \
                  --cpus {resources.threads} \
                  1> {log.o} \
                  2> {log.e}
                """

rule kraken2:
    input:
                "raw_data/fna/{isolate}.fna",
    output:
                "analysis/kraken2/kraken2_{isolate}.out",
    envmodules:
                "kraken2/2.1.3",
    log:
                "logs/kraken2/kraken2_{isolate}.log",
    resources:
                threads = 12,
                time = "12:00:00",
                nodes = 1,
                mem_gb = 252,
                partition = "pgec-mem",
                qos = "pgec",
                name = "kraken2",
    shell:
                """
                kraken2 \
                    --db /reference/data/NCBI/kraken2/2023-09-10/pluspfp/ \
                    --threads {resources.threads} \
                    --output {output} \
                    --use-names \
                    {input} \
                    2>&1 {log}
                """

rule kraken2_concat:
    input:
                "raw_data/fna_concat/{isolate}_concat.fna",
    output:
                "analysis/kraken2/kraken2_concat_{isolate}.out",
    envmodules:
                "kraken2/2.1.3",
    log:
                "logs/kraken2/kraken2_concat_{isolate}.log",
    resources:
                threads = 12,
                time = "12:00:00",
                nodes = 1,
                mem_gb = 252,
                partition = "pgec-mem",
                qos = "pgec",
                name = "kraken2",
    shell:
                """
                kraken2 \
                    --db /reference/data/NCBI/kraken2/2023-09-10/pluspfp/ \
                    --threads {resources.threads} \
                    --output {output} \
                    --use-names \
                    {input} \
                    2>&1 {log}
                """

rule kraken2_summarize:
    input:
                expand("analysis/kraken2/kraken2_{isolate}.out", isolate=meta["isolate"].tolist())
    output:
                "analysis/kraken2/kraken2_summarize.txt",
    log:
                "logs/kraken2/kraken2_summarize.log",
    resources:
                threads = 1,
                time = "00:10:00",
                nodes = 1,
                mem_gb = 8,
                partition = "pgec-mem",
                qos = "pgec",
                name = "kraken2_summarize",
    shell:
                """
                awk '{{print $1,$2,$3,$4,$5}}' analysis/kraken2/kraken2_*.out > {output}
                """

rule kraken2_concat_summarize:
    input:
                expand("analysis/kraken2/kraken2_concat_{isolate}.out", isolate=meta["isolate"].tolist())
    output:
                "analysis/kraken2/kraken2_concat_summarize.txt",
    log:
                "logs/kraken2/kraken2_concat_summarize.log",
    resources:
                threads = 11,
                time = "00:10:00",
                nodes = 1,
                mem_gb = 8,
                partition = "pgec-mem",
                qos = "pgec",
                name = "kraken2_summarize",
    shell:
                """
                awk '{{print $1,$2,$3,$4,$5}}' analysis/kraken2/kraken2_concat_*.out > {output}
                """