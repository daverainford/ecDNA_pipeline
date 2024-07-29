t_bam = list.files(path = "./", pattern = "_T\\.bam$")
n_bam = list.files(path = ".", pattern = "_N\\.bam")
split = strsplit(t_bam, "_")
sample_id = sapply(split, function(x) x[1])

gcap.workflow(
  tumourseqfile = t_bam, normalseqfile = n_bam, 
  tumourname = paste0(sample_id, "_T"), normalname = paste0(sample_id, "_N"), jobname = sample_id,
  outdir = "./",
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter", 
  g1000allelesprefix = file.path(
    "reference/",
    "G1000_allelesAll_hg38"
  ), 
  g1000lociprefix = file.path("reference/",
                              "G1000_lociAll_hg38"
  ),
  GCcontentfile = "reference/GC_G1000_hg38.txt",
  replictimingfile = "reference/RT_G1000_hg38.txt",
  skip_finished_ASCAT = TRUE,
  skip_ascat_call = FALSE,
  result_file_prefix = sample_id,
  extra_info = df,
  include_type = FALSE,
  genome_build = "hg38",
  model = "XGB11"
)
