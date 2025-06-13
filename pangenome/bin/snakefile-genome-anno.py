import pandas as pd
import os
from glob import glob
from shutil import which
from snakemake.utils import min_version

### set minimum snakemake version ###
min_version("8.25.0")

meta = pd.read_csv("bin/genome_metadata.tsv", sep="\t")
# meta["isolate"].tolist()

Ns = [0,10,20,30,40,50,60,70,80,90,100]

rule all:
    input:
                # gtdbtk
                "analysis/gtdbtk/",
                # fastani
                "analysis/fastani/fastani_many-to-many.tsv",
                # gtotree
                "analysis/gtotree/",
                # gtotree w/ Nocardiodes
                "analysis/gtotreeN/",
                # prokka
                expand("analysis/prokka/{isolate}.{ext}", isolate=meta["isolate"].tolist(), ext = ["gff","gbk","fna","faa","ffn","sqn","fsa","tbl","err","log","txt","tsv"]),
                # orthofinder
                "analysis/orthofinder/",
                # quast_N
                expand("analysis/quast_N/quast_N{N}/report.tsv", N = [0,10,20,30,40,50,60,70,80,90,100]),
                # get_Nx
                "analysis/quast_Nx.tsv",
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

rule fastani:
    input:
                "bin/fastani-genome-list.txt",
    output:
                tsv = "analysis/fastani/fastani_many-to-many.tsv",
                mat = "analysis/fastani/fastani_many-to-many.tsv.matrix",
    conda:
                "envs/fastani.yaml",
    resources:
                threads = 8,
                time = "48:00:00",
                nodes = 1,
                mem_gb = 32,
                partition = "short",
                qos = "pgec",
                name = "fastani",
    log:
                out = "logs/fastani.o",
                err = "logs/fastani.e",
    shell:
                """
                fastANI \
                    --ql {input} \
                    --rl {input} \
                    --threads {resources.threads} \
                    --matrix \
                    -o {output.tsv} \
                    2> {log.out} \
                    1> {log.err}
                """

rule gtotree:
    input:
                paths = "bin/gtotree-genomes.tsv",
                labels = "bin/gtotree-labels.tsv",
    output:
                directory("analysis/gtotree/"),
    conda:
                "envs/gtotree.yaml",
    envmodules:
                "python_3/3.11.1",
    resources:
                threads = 48,
                time = "48:00:00",
                nodes = 1,
                mem_gb = 120,
                partition = "priority-mem",
                qos = "pgec",
                name = "gtotree",
    log:
                out = "logs/gtotree.o",
                err = "logs/gtotree.e",
    shell:
                """
                GToTree \
                    -f {input.paths} \
                    -m {input.labels} \
                    -o {output} \
                    -H Actinobacteria \
                    -j {resources.threads} \
                    1> {log.out} \
                    2> {log.err}
                """

rule gtotreeN:
    input:
                paths = "bin/gtotreeN-genomes.tsv",
                labels = "bin/gtotreeN-labels.tsv",
    output:
                directory("analysis/gtotreeN/"),
    conda:
                "envs/gtotree.yaml",
    envmodules:
                "python_3/3.11.1",
    resources:
                threads = 48,
                time = "48:00:00",
                nodes = 1,
                mem_gb = 120,
                partition = "priority-mem",
                qos = "pgec",
                name = "gtotree",
    log:
                out = "logs/gtotree.o",
                err = "logs/gtotree.e",
    shell:
                """
                GToTree \
                    -f {input.paths} \
                    -m {input.labels} \
                    -o {output} \
                    -H Actinobacteria \
                    -j {resources.threads} \
                    1> {log.out} \
                    2> {log.err}
                """

rule prokka:
    input:
                "raw_data/fna/{isolate}.fna",
    output:
                "analysis/prokka/{isolate}.gff",
                "analysis/prokka/{isolate}.gbk",
                "analysis/prokka/{isolate}.fna",
                "analysis/prokka/{isolate}.faa",
                "analysis/prokka/{isolate}.ffn",
                "analysis/prokka/{isolate}.sqn",
                "analysis/prokka/{isolate}.fsa",
                "analysis/prokka/{isolate}.tbl",
                "analysis/prokka/{isolate}.err",
                "analysis/prokka/{isolate}.log",
                "analysis/prokka/{isolate}.txt",
                "analysis/prokka/{isolate}.tsv",
    conda:
                "envs/prokka.yaml",
    resources:
                threads = 8,
                time = "48:00:00",
                nodes = 1,
                mem_gb = 32,
                partition = "short",
                qos = "pgec",
                name = "prokka",
    log:
                out = "logs/prokka/prokka_{isolate}.o",
                err = "logs/prokka/prokka_{isolate}.e",
    shell:
                """
                prokka \
                --outdir analysis/prokka/ \
                --force \
                --prefix {wildcards.isolate} \
                {input} \
                1> {log.out} \
                2> {log.err}
                """

rule orthofinder:
    input:
                expand("raw_data/faa/{isolate}.faa", isolate = meta["isolate"].tolist()),
    params:
                in_dir = "raw_data/faa/",
    output:
                directory("analysis/orthofinder/results_orthofinder"),
                # analysis/orthofinder/results_orthofinder/
    conda:
                "envs/orthofinder.yaml",
    resources:
                threads = 48,
                time = "48:00:00",
                nodes = 1,
                mem_gb = 128,
                partition = "priority-mem",
                qos = "pgec",
                name = "orthofinder",
    log:
                out = "logs/orthofinder.o",
                err = "logs/orthofinder.e",
    shell:
                """
                orthofinder \
                    -a {resources.threads} \
                    -f {params.in_dir} \
                    -T fasttree \
                    -n orthofinder \
                    -o {output} \
                    1> {log.out} \
                    2> {log.err}
                """

rule quast_N:
    input:
                expand("raw_data/fna/{isolate}.fna", isolate = meta["isolate"].tolist()),
    params:
                out_dir = "analysis/quast_N/quast_N{N}",
    output:
                "analysis/quast_N/quast_N{N}/report.tsv",
    conda:
                "envs/quast.yaml",
    resources:
                threads = 8,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 16,
                partition = "short",
                qos = "pgec",
                name = "quast",
    log:
                "logs/quast_N/quast_N{N}",
    shell:
                """
                quast \
                    -o {params.out_dir} \
                    --x-for-Nx {wildcards.N} \
                    {input}
                """

rule get_Nx:
    input:
                expand("analysis/quast_N/quast_N{N}/report.tsv", N =[0,10,20,30,40,50,60,70,80,90,100]),
    resources:
                threads = 8,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 16,
                partition = "short",
                qos = "pgec",
                name = "quast",
    output:
                "analysis/quast_Nx.tsv",
    shell:
                """
                Rscript \
                    --vanilla \
                    bin/get_Nx.R
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
