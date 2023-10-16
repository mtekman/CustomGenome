sub4 = system.file("extdata", "subsample4.gtf.gz", package="CustomGenome")
sub8 = system.file("extdata", "subsample8.gtf.gz", package="CustomGenome")
usfa = system.file("extdata", "user_sequences.fa", package="CustomGenome")

test_get_genome_urls_default <- function() {
  checkEquals(get_genome_urls(),
              list(gtf = paste0("http://ftp.ensembl.org/pub/release-105/gtf/",
                                "mus_musculus/Mus_musculus.GRCm39.105.gtf.gz"),
                   fasta = paste0("http://ftp.ensembl.org/pub/release-105/fasta/",
                                  "mus_musculus/dna/Mus_musculus.GRCm39.",
                                  "dna.primary_assembly.fa.gz")))
}

test_get_genome_urls_danrer_release104 <- function() {
  checkEquals(get_genome_urls(species="danio_rerio", release="104"),
              list(gtf = paste0("http://ftp.ensembl.org/pub/release-104/gtf/",
                                "danio_rerio/Danio_rerio.GRCz11.104.gtf.gz"),
                   fasta = paste0("http://ftp.ensembl.org/pub/release-104/fasta/",
                                  "danio_rerio/dna/Danio_rerio.GRCz11.",
                                  "dna.primary_assembly.fa.gz")))
}

test_get_genome_urls_custom <- function() {
  checkEquals(get_genome_urls(release="80", build="EquCab2", species="equus_caballus",
                              fasta_type="dna.toplevel"),
              list(gtf = paste0("http://ftp.ensembl.org/pub/release-80/gtf/",
                                "equus_caballus/Equus_caballus.EquCab2.80.gtf.gz"),
                   fasta = paste0("http://ftp.ensembl.org/pub/release-80/fasta/",
                                  "equus_caballus/dna/Equus_caballus.EquCab2.",
                                  "dna.toplevel.fa.gz")))
}

test_get_genome_files <- function() {
  checkEquals(get_genome_files(
    fasta_type = "dna_rm.nonchromosomal",
    gtf_type = "abinitio.gtf", output_folder="/tmp"),
    list(gtf = "/tmp/Mus_musculus.GRCm39.105.abinitio.gtf.gz",
         fasta = "/tmp/Mus_musculus.GRCm39.dna_rm.nonchromosomal.fa.gz"))
}

## Dependent on last test running
test_check_sum_matches1 <- function() {
  checkEquals(check_sum_matches(
    paste0("http://ftp.ensembl.org/pub/release-105/gtf/",
           "mus_musculus/Mus_musculus.GRCm39.105.abinitio.gtf.gz"),
    "/tmp/Mus_musculus.GRCm39.105.abinitio.gtf.gz"), TRUE)
}

test_check_sum_matches2 <- function() {
  checkEquals(check_sum_matches(
    paste0("http://ftp.ensembl.org/pub/release-105/fasta/",
           "mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.nonchromosomal.fa.gz"),
    "/tmp/Mus_musculus.GRCm39.dna_rm.nonchromosomal.fa.gz"), TRUE)
}

test_qc_filter_lines_of_gtf_1 <- function() {
  checkEquals(sapply(
    strsplit(gsub("(\n| )", "",
                  capture_messages(
                    qc_filter_lines_of_gtf(sub8)))[2:4], split=":"),
    function(x) as.integer(x[[2]])),
    c(100, 80, 78))
}

test_qc_filter_lines_of_gtf_2 <- function() {
  checkEquals(sapply(
    strsplit(gsub("(\n| )", "",
                  capture_messages(
                    qc_filter_lines_of_gtf(sub4)))[2:4], split=":"),
    function(x) as.integer(x[[2]]))
    c(100, 88, 86))
}

test_add_lines_to_gtf <- function() {
  obj_gtf = as.data.frame(import(sub8))
  obj_fas = readDNAStringSet(usfa)
  add_seqs_to_gtf(obj_gtf, obj_fas)
  checkEquals(TRUE, FALSE)
}

test_add_seqs_to_gtf_and_fasta <- function() {
  checkEquals(TRUE, FALSE)
}

test_prime_fastq_files <- function() {
  checkEquals(TRUE, FALSE)
}

test_retrieve_index <- function() {
  checkEquals(TRUE, FALSE)
}

test_perform_alignment <- function() {
  checkEquals(TRUE, FALSE)
}

test_summarize_alignment <- function() {
  checkEquals(TRUE, FALSE)
}

test_generate_count_matrix <- function() {
  checkEquals(TRUE, FALSE)
}
