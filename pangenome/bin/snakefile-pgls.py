import pandas as pd
import os
from glob import glob
from shutil import which
from snakemake.utils import min_version

### set minimum snakemake version ###
min_version("8.25.0")

# change this to test additional phenotype associations
PHENES = ["de"]

SPLIT_N = list(range(1,28,1))

rule all:
    input:
                # split
                "analysis/pgls_split/",
                # pgls
                # expand("analysis/pgls2/{phene}/{phene}-{split_n}.tsv", phene = PHENES, split_n = SPLIT_N),
                # pgls_merge
                # expand("analysis/pgls2/results/{phene}.tsv", phene = PHENES),
                # mafft
                # expand("analysis/mafft/{OG}.fa", OG = "OG0000000"),
                # hmmbuild_hmmemit
                # "analysis/hmmbuild/{OG}.hmm",
                # "analysis/hmmbuild/{OG}.txt",
                # "analysis/hmmemit/{OG}.fa",
                # hmmemit_consensus
                "analysis/hmmemit_consensus/consensus.fa"
                # pgls_annotate
                expand("analysis/pgls/pgls_annotate/{phene}-anno.tsv", phene = PHENES),

checkpoint split:
    input:
                og = "analysis/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv",
                col = "raw_data/colonization.tsv",
                plant = "raw_data/microbox-plant-phenes_mean.tsv",
    params:
                script = "bin/split.R",
    output:
                out_dir = directory("analysis/pgls_split/"),
    conda:
                "envs/r.yaml",
    log:
                o = "logs/split.o",
                e = "logs/split.e"
    resources:
                threads = 1,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 8,
                partition = "priority-mem",
                qos = "pgec",
                name = "split",
    shell:
                """
                Rscript \
                  --vanilla \
                  {params.script} \
                  {input.og} \
                  {input.col} \
                  {input.plant} \
                  {output.out_dir} \
                  1> {log.o} \
                  2> {log.e}
                """

rule pgls2:
    input:
                split_file = "analysis/pgls_split/split-{split_n}.tsv",
                tree_file =  "analysis/gtotreeN/gtotreeN.tre",
    params:
                script = "bin/pgls2.R",
                phene =  "{phene}",
    output:
                "analysis/pgls2/{phene}/{phene}-{split_n}.tsv",
    conda:
                "envs/r.yaml",
    log:
                "logs/pgls2/{phene}/{phene}-{split_n}.o",
    resources:
                threads = 1,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 16,
                partition = "priority-mem",
                qos = "pgec",
                name = "pgls2",
    shell:
                """
                Rscript \
                  --vanilla \
                  {params.script} \
                  {input.split_file} \
                  {input.tree_file} \
                  {params.phene} \
                  {output} \
                  1> {log} \
                  2> {log}
                """

def aggregate_pgls(wildcards):
    checkpoint_output = checkpoints.split.get(**wildcards).output[0]
    x = expand("analysis/pgls2/{phene}/{phene}-{split_n}.tsv",
        phene=wildcards.phene,
        split_n=glob_wildcards(os.path.join(checkpoint_output, "split-{i}.tsv")).i)
    return x

checkpoint pgls_merge:
    input:
                aggregate_pgls,
    params:
                script = "bin/pgls_merge.R",
                phene = "{phene}",
                in_dir = "analysis/pgls2/{phene}",
                cutoff = 0.25,
    output:
                "analysis/pgls2/results/{phene}.tsv",
    conda:
                "envs/r.yaml",
    log:
                "logs/pgls_merge/pgls_merge-{phene}.o",
    resources:
                threads = 1,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 16,
                partition = "priority-mem",
                qos = "pgec",
                name = "pgls_merge",
    shell:
                """
                Rscript \
                    --vanilla \
                    {params.script} \
                    {params.phene} \
                    {params.in_dir} \
                    {params.cutoff} \
                    {output} \
                    1> {log} \
                    2> {log}
                """

rule mafft:
  input:
                "analysis/orthofinder/Orthogroup_Sequences/{OG}.fa",
  output:
                "analysis/mafft/{OG}.fa",
  log:
                "logs/mafft/mafft-{OG}.log",
  envmodules:
                "mafft/7.505",
  resources:
                threads = 8,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 16,
                partition = "priority-mem",
                qos = "pgec",
                name = "mafft",
  shell:
                """
                mafft \
                  --thread {resources.threads} \
                  --auto \
                  {input} \
                  > {output} \
                  2> {log}
                """

rule hmmbuild_hmmemit:
  input:
                # pgls = "analysis/pgls/results/de.tsv",
                msa = "analysis/mafft/{OG}.fa",
  output:
                hmm = "analysis/hmmbuild/{OG}.hmm",
                txt = "analysis/hmmbuild/{OG}.txt",
                cons = "analysis/hmmemit/{OG}.fa",
  log:
                "logs/hmmbuild/hmmbuild-{OG}.o"
  envmodules:
                "hmmer3/3.3.2",
  resources:
                threads = 4,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 8,
                partition = "priority-mem",
                qos = "pgec",
                name = "hmmbuild_emit",
  shell:
                """
                hmmbuild \
                  --amino \
                  -n {wildcards.OG} \
                  -o {output.txt} \
                  {output.hmm} \
                  {input.msa} \
                  1> {log} \
                  2> {log}

                hmmemit \
                  -c \
                  {output.hmm} \
                  1> {output.cons} \
                  2> {log}
                """

def hmmemit_join(wildcards):
    df = pd.read_csv("analysis/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv", sep="\t")
    OGs = df["orthogroup"].tolist()
    files = expand("analysis/hmmemit/{OG}.fa", OG = OGs)
    return files

rule hmmemit_consensus:
  input:
                hmmemit_join
  output:
                "analysis/hmmemit_consensus/consensus.fa",
  resources:
                threads = 4,
                time = "1:00:00",
                nodes = 1,
                mem_gb = 8,
                partition = "priority-mem",
                qos = "pgec",
                name = "consensus",
  shell:
                """
                cat {input} >> {output}
                """

rule pgls_annotate:
  input:
              pgls = "analysis/pgls/results/{phene}.tsv",
              egg = "analysis/emapper/out.emapper.annotations",
              fun = "symlink_data/cog-2020/fun-20.tab",
  params:
              script = "bin/pgls_annotate.R",
  output:
              "analysis/pgls/pgls_annotate/{phene}-anno.tsv",
  log:
              "logs/pgls_annotate/pgls_annotate-{phene}.log",
  shell:
              """
              Rscript \
                --vanilla \
                {params.script} \
                {input.pgls} \
                {input.egg} \
                {input.fun} \
                {output} \
                1> {log} \
                2> {log}
              """

## DO NOT INCLUDE

# def hmmer_result(wildcards):
#     df = pd.read_csv("analysis/pgls/results/de.tsv", sep="\t")
#     OGs = df["orthogroup"].tolist()
#     files = expand("analysis/hmmer/{OG}-{db}.tbl", OG = OGs, db = ["pfam","uniref"])
#     return files

# rule hmmer_merge:
#   input:
#               # hmmer_result,
#               "analysis/hmmer/OG0000006-uniref.tbl",
#   params:
#               script = "bin/hmmer_merge.R",
#   output:
#               "analysis/hmmer_merge.tsv",
#   log:
#               "logs/hmmer_merge.log",
#   shell:
#               """
#               Rscript \
#                 --vanilla \
#                 {params.script} \
#                 {output} \
#                 1> {log} \
#                 2> {log}
#               """