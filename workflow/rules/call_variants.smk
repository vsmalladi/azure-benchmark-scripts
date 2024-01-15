#!/usr/bin/env snakemake

"""
Call variants relative to a reference genome
"""

#### Output ####
HAPLOTYPER = "haplotyper/{sample}/calls.vcf.gz"
DNASCOPE_TMP = "dnascope_tmp/{sample}/calls.vcf.gz"
DNASCOPE = "dnascope/{sample}/calls.vcf.gz"
DNASCOPE_LR = "dnascope_lr/{sample}/calls.vcf.gz"

#### Functions ####
# Get the read filter for Ultima data
def get_dnascope_rf(wldc, input):
  if "Ultima" not in wldc.sample:
    return ""
  return f"--read_filter UltimaReadFilter"

# Get the proper alignemnt file for calling
def get_alignemnts(get_idx=False):
  def _get_alignemnts(wldc):
    if "Ultima" in wldc.sample:
      fn = CRAM.format(sample=wldc.sample)
      return fn + ".crai" if get_idx else fn
    else:
      fn = DEDUP.format(sample=wldc.sample, bwa_model="yes")
      return fn + ".bai" if get_idx else fn
  return _get_alignemnts

def get_tech(wldc):
  if "PacBio" in wldc.sample:
    return "HiFi"
  else:
    return "ONT"

def get_dnascopelr_bams(suffix=""):
  def _get_dnascopelr_bams(wldc):
    expected = []
    sample_list = config["input"]["samples"][wldc.sample]

    n_splits = 1
    with checkpoints.numa_splits.get().output[0].open() as fh:
      n_splits = int(fh.read().rstrip())

    for read_idx in range(len(sample_list)):
      for subset in range(n_splits):
        expected.append(MINIMAP2.format(subset=subset, n_subsets=n_splits, read_idx=read_idx, sample=wldc.sample) + suffix)
    return expected
  return _get_dnascopelr_bams


#### Rules ####
rule haplotyper:
  input:
    bam = lambda wldc: DEDUP.format(sample=wldc.sample, bwa_model="none"),
    bai = lambda wldc: DEDUP.format(sample=wldc.sample, bwa_model="none") + ".bai",
    qcal = rules.qualcal.output,
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    autosomes_bed = rules.autosomes_bed.output,
    dbsnp = SITE_VCF.format(vcf_name="dbsnp"),
    sentieon = config["tools"]["sentieon"],
  output:
    vcf = HAPLOTYPER,
    tbi = HAPLOTYPER + ".tbi",
  benchmark:
    HAPLOTYPER + ".benchmark.txt"
  log:
    stdout = HAPLOTYPER + ".stdout",
    stderr = HAPLOTYPER + ".stderr",
  params:
    xargs = "--pcr_indel_model none"
  threads:
    300
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver -i "{input.bam}" -r "{input.fa}" -q "{input.qcal}" \
      --interval "{input.autosomes_bed}" -t {threads} --algo Haplotyper -d "{input.dbsnp}" {params.xargs} "{output.vcf}"
    """

rule dnascope:
  input:
    bam = get_alignemnts(),
    bai = get_alignemnts(True),
    fa = get_qualcal_ref(),
    fai = get_qualcal_ref(".fai"),
    autosomes_bed = rules.autosomes_bed.output,
    model = get_dnascope_model,
    dbsnp = SITE_VCF.format(vcf_name="dbsnp"),
    sentieon = config["tools"]["sentieon"],
  output:
    vcf = DNASCOPE_TMP,
    tbi = DNASCOPE_TMP + ".tbi",
  benchmark:
    DNASCOPE_TMP + ".benchmark.txt"
  log:
    stdout = DNASCOPE_TMP + ".stdout",
    stderr = DNASCOPE_TMP + ".stderr",
  params:
    xargs = "--pcr_indel_model none",
    read_filter = get_dnascope_rf,
  threads:
    300
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver -i "{input.bam}" -r "{input.fa}" {params.read_filter} \
      --interval "{input.autosomes_bed}" -t {threads} --algo DNAscope {params.xargs} \
      --model "{input.model}/dnascope.model" -d "{input.dbsnp}" "{output.vcf}"
    """

rule dnamodelapply:
  input:
    vcf = rules.dnascope.output.vcf,
    tbi = rules.dnascope.output.tbi,
    fa = get_qualcal_ref(),
    fai = get_qualcal_ref(".fai"),
    model = get_dnascope_model,
    sentieon = config["tools"]["sentieon"],
  output:
    vcf = DNASCOPE,
    tbi = DNASCOPE + ".tbi",
  benchmark:
    DNASCOPE + ".benchmark.txt"
  log:
    stdout = DNASCOPE + ".stdout",
    stderr = DNASCOPE + ".stderr",
  threads:
    300
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    {input.sentieon} driver -t {threads} -r "{input.fa}" --algo DNAModelApply \
      --model "{input.model}/dnascope.model" -v "{input.vcf}" "{output.vcf}"
    """

rule dnascope_longread:
  input:
    bams = get_dnascopelr_bams(),
    bais = get_dnascopelr_bams(".bai"),
    fa = rules.download_ref.output.fa,
    fai = rules.download_ref.output.fai,
    autosomes_bed = rules.autosomes_bed.output,
    model = get_dnascope_model,
    dbsnp = SITE_VCF.format(vcf_name="dbsnp"),
    sentieon_cli = config["tools"]["sentieon_cli"],
    sentieon = config["tools"]["sentieon"],
    bedtools = config["tools"]["bedtools"],
    bcftools = config["tools"]["bcftools"],
  output:
    vcf = DNASCOPE_LR,
    tbi = DNASCOPE_LR + ".tbi",
  benchmark:
    DNASCOPE_LR + ".benchmark.txt"
  log:
    stdout = DNASCOPE_LR + ".stdout",
    stderr = DNASCOPE_LR + ".stderr",
  params:
    tech = get_tech,
    input = lambda wldc, input: " '" + "' '".join(input.bams) + "'",
  threads:
    300
  shell:
    """
    set -exvuo pipefail
    outdir=$(dirname "{output.vcf}")
    mkdir -p "$outdir"
    exec 1>"{log.stdout}" 2>"{log.stderr}"

    sentieon_dir=$(dirname "{input.sentieon}")
    bedtools_dir=$(dirname "{input.bedtools}")
    bcftools_dir=$(dirname "{input.bcftools}")

    export PATH=$sentieon_dir:$bedtools_dir:$bcftools_dir:$PATH
    "{input.sentieon_cli}" dnascope-longread --tech "{params.tech}" \
      -r "{input.fa}" -i {params.input} -m "{input.model}" \
      -b "{input.autosomes_bed}" -d "{input.dbsnp}" "{output.vcf}"
    """
