## Debug
## sub4 = file.path("inst/extdata", "subsample4.gtf.gz")
## sub8 = file.path("inst/extdata", "subsample8.gtf.gz")
## usfa = file.path("inst/extdata", "user_sequences.fa")
## tiny = file.path("inst/extdata", "tiny.fa.gz")

sub4 <- system.file("extdata", "subsample4.gtf.gz", package="CustomGenome")
sub8 <- system.file("extdata", "subsample8.gtf.gz", package="CustomGenome")
usfa <- system.file("extdata", "user_sequences.fa", package="CustomGenome")
tiny <- system.file("extdata", "tiny.fa.gz", package="CustomGenome")
tiny_gtf <- system.file("extdata", "tiny.gtf.gz", package="CustomGenome")
r1_fastq <- system.file("extdata", "r1.fq.gz", package="CustomGenome")
r2_fastq <- system.file("extdata", "r2.fq.gz", package="CustomGenome")

test_that("sysfiles", {
  expect_equal(file.exists(sub4), TRUE)
  expect_equal(file.exists(sub8), TRUE)
  expect_equal(file.exists(usfa), TRUE)
  expect_equal(file.exists(tiny), TRUE)
})

tmp_dir <- tempdir()

test_that("get_genome_urls_default", {
  expect_equal(
    get_genome_urls(),
    list(gtf = paste0("http://ftp.ensembl.org/pub/release-105/gtf/",
                      "mus_musculus/Mus_musculus.GRCm39.105.gtf.gz"),
         fasta = paste0("http://ftp.ensembl.org/pub/release-105/fasta/",
                        "mus_musculus/dna/Mus_musculus.GRCm39.",
                        "dna.primary_assembly.fa.gz")))
})

test_that("get_genome_urls_danrer_release104", {
  expect_equal(
    get_genome_urls(species="danio_rerio", release="104"),
    list(gtf = paste0("http://ftp.ensembl.org/pub/release-104/gtf/",
                      "danio_rerio/Danio_rerio.GRCz11.104.gtf.gz"),
         fasta = paste0("http://ftp.ensembl.org/pub/release-104/fasta/",
                        "danio_rerio/dna/Danio_rerio.GRCz11.",
                        "dna.primary_assembly.fa.gz")))
})

test_that("get_genome_urls_custom", {
  expect_equal(
    get_genome_urls(release="80", build="EquCab2", species="equus_caballus",
                    fasta_type="dna.toplevel"),
    list(gtf = paste0("http://ftp.ensembl.org/pub/release-80/gtf/",
                      "equus_caballus/Equus_caballus.EquCab2.80.gtf.gz"),
         fasta = paste0("http://ftp.ensembl.org/pub/release-80/fasta/",
                        "equus_caballus/dna/Equus_caballus.EquCab2.",
                        "dna.toplevel.fa.gz")))
})

test_that("get_genome_files", {
  capture_messages(res <- get_genome_files(
                     fasta_type = "dna_rm.nonchromosomal",
                     gtf_type = "abinitio.gtf", output_folder = tmp_dir
                   ))
  expect_equal(res,
               list(
                 gtf = file.path(tmp_dir, "Mus_musculus.GRCm39.105.abinitio.gtf.gz"),
                 fasta = file.path(tmp_dir,
                                   "Mus_musculus.GRCm39.dna_rm.nonchromosomal.fa.gz")
               ))
})

## Dependent on last test running
test_that("check_sum_matches", {
  capture_messages(res <- check_sum_matches(
                     paste0(
                       "http://ftp.ensembl.org/pub/release-105/gtf/",
                       "mus_musculus/Mus_musculus.GRCm39.105.abinitio.gtf.gz"
                     ),
                     file.path(tmp_dir, "Mus_musculus.GRCm39.105.abinitio.gtf.gz")
                   ))
  expect_equal(res, TRUE)
})

test_that("check_sum_matches", {
  capture_messages(
    res <- check_sum_matches(
      paste0(
        "http://ftp.ensembl.org/pub/release-105/fasta/",
        "mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.nonchromosomal.fa.gz"
      ),
      file.path(tmp_dir, "Mus_musculus.GRCm39.dna_rm.nonchromosomal.fa.gz")
    ))
  expect_equal(res, TRUE)
})

test_that("qc_filter_lines_of_gtf_1", {
  expect_equal(sapply(
    strsplit(gsub("(\n| )", "",
                  capture_messages(
                    qc_filter_lines_of_gtf(sub8)))[2:4], split=":"),
    function(x) as.integer(x[[2]])), c(100, 80, 78))
})

test_that("qc_filter_lines_of_gtf_2", {
  expect_equal(sapply(
    strsplit(gsub("(\n| )", "",
                  capture_messages(
                    qc_filter_lines_of_gtf(sub4)))[2:4], split=":"),
    function(x) as.integer(x[[2]])), c(100, 88, 86))
})

test_that("add_lines_to_gtf", {
  obj_gtf <- import(sub8)
  obj_fas <- readDNAStringSet(usfa)
  capture_messages(res <- add_lines_to_gtf(obj_gtf, obj_fas))
  expect_equal(
    as.character(tail(res, 3)@seqnames@values),
    c("One", "Two", "Three")
  )
})

## This is purely used to generate test-data
generate_tiny_genome <- function(genome_file,
                                 n_samp = 10, cont_stretch = 1e4) {
  gen <- readDNAStringSet(genome_file)
  gen <- head(gen, 21)
  for (chrom in names(gen)) {
    message("Processing ", chrom)
    max_cap = (length(gen[[chrom]]) - (cont_stretch + 10))
    gen[[chrom]] <- gen[[chrom]][sort(
                         do.call(c, lapply(
                                      sample(1:max_cap, n_samp),
                                      function(x) {
                                        seq(x, x + cont_stretch - 1)
                                      })))]
  }
  names(gen) <- gsub("^([0-9A-Z]+)\\s+.*", "\\1", names(gen))
  writeXStringSet(gen, "tiny.fa.gz", compress = TRUE)
  return(gen)
}

test_that("add_seqs_to_gtf_and_fasta", {
  gfiles <- list(gtf = sub8, fasta = tiny)
  capture_messages(res <- add_seqs_to_gtf_and_fasta(gfiles, usfa))
  expect_equal(
    as.character(tail(import(res$gtf), 3)@seqnames@values),
    c("One", "Two", "Three")
  )
  expect_equal(
    tail(names(readDNAStringSet(res$fasta)), 3),
    c("One", "Two", "Three")
  )
  file.remove(unlist(res)) ## remove to not interfere with other tests
})

test_that("prime_fastq_files_1", {
  ## Paired
  expect_equal(
    prime_fastq_files(
      system.file("extdata", package = "CustomGenome"),
      "e4.gtf.gz", "e8.gtf.gz", "e-align.bam"
    ),
    list(reads1 = sub4, reads2 = sub8, align_base = "subsample-align.bam")
  )
  ## Single
  mes <- capture_messages(res <- prime_fastq_files(
                            system.file("extdata",
                                        package = "CustomGenome"), ".gtf.gz",
                            NULL, ".bam"))
  expect_equal(
    res,
    list(reads1 = c(sub4, sub8, tiny_gtf), reads2 = NULL,
         align_base = sub(".gtf.gz", ".bam", basename(c(sub4, sub8, tiny_gtf))))
  )
  expect_equal(mes, "We assume that these are single-end reads\n")
})

test_that("prime_fastq_files_2", {
  tdir <- tempdir()
  fnames <- file.path(tdir,
                      c("sample1_r1.fa.gz", "sample1_r2.fa.gz",
                        "sample1_r3.fa.gz",
                        "s2_r1.fa.gz", "s2_r2.fa.gz", "s2_i3.fa.gz",
                        "reads_r1.fa.gz", "reads_r2.fa.gz"))
  file.create(fnames)

  expect_equal(
    prime_fastq_files(tdir, "_r1.fa.gz", "_r2.fa.gz", "-align.bam"),
    list(
      reads1 = file.path(tdir, c("reads_r1.fa.gz", "s2_r1.fa.gz",
                                 "sample1_r1.fa.gz")),
      reads2 = file.path(tdir, c("reads_r2.fa.gz", "s2_r2.fa.gz",
                                 "sample1_r2.fa.gz")),
      align_base = c("reads-align.bam", "s2-align.bam", "sample1-align.bam")
    )
  )
  file.remove(fnames)
})

test_that("all_subread_commands", {
  require(Rsamtools)
  test_dir <- tempdir()

  read_lists <- list(reads1 = c(r1_fastq), reads2 = c(r2_fastq),
                     align_base = "test.bam")
  dir_lists <- list(
    index = paste0(test_dir, "-subread-index"),
    align = paste0(test_dir, "-aligned"),
    count = paste0(test_dir, "-counts"),
    stats = paste0(test_dir, "-stats")
  )

  if (dir.exists(dir_lists$index)) {
    unlink(dir_lists$index, recursive = TRUE)
  }
  ## Build Index
  mes <- capture_messages(capture_output(
    zz <- retrieve_index(tiny, dir_lists$index)))
  ## Check Index
  expect_equal(
    readLines(file.path(dir_lists$index, "reference_index.reads"), n = 5)[[5]],
    "512020\t13"
  )
  expect_equal(substr(mes, 1, 14), "Building Index")

  ## Run again, this time with retrieval only
  mes <- capture_messages(res <- retrieve_index(tiny, dir_lists$index))
  expect_equal(res, file.path(dir_lists$index, "reference_index"))
  expect_equal(mes, paste0("Index already built: ", dir_lists$index, "\n"))

  ## Perform Alignment
  dir_lists$index <- res
  mes <- capture_output(res <- perform_alignment(dir_lists, read_lists,
                                                 nthreads = 1))
  ## Check Alignment
  expect_equal(res["Mapped_fragments", ], 1)
  expect_equal(gsub("\\s\\s+", "", unlist(strsplit(mes, split = "\n"))[62]),
               "||Not properly paired : 1||")

  ## Check Alignment File
  bamfile <- file.path(dir_lists$align, read_lists$align_base)
  expect_equal(
    as.character(scanBam(bamfile)[[1]]$seq[[1]][1:10]),
    "ATATTAGTTG")

  ## Check summary
  mes <- capture_messages(res <- summarize_alignment(
                            dir_lists,
                            read_lists$align_base))
  expect_equal(
    file.exists(file.path(dir_lists$align, "test.bam.summary")),
    TRUE
  )
  expect_equal(res$Percentage_Mapped, 100)

  ## Do quantification
  cfiles <- generate_count_matrix(dir_lists,
                                  gtf_file = tiny_gtf,
                                  bam_pattern = "*test.bam$",
                                  GTF.attrType = "gene_id", nthreads = 30)
  ## Test quantification
  expect_equal(readLines(cfiles$stat)[2], "Assigned\t1")
  expect_equal(readLines(cfiles$count)[3], "large_hit\t1")

  ## Clean up
  if (dir.exists(dir_lists$index)) {
    unlink(dir_lists$index, recursive = TRUE)
  }
  if (dir.exists(dir_lists$align)) {
    unlink(dir_lists$index, recursive = TRUE)
  }
  unlink(tempdir(), recursive=TRUE)
})
