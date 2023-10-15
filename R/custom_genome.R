####################################################################
## required packages
####################################################################
require(rtracklayer) ## For reading GTF via `import'
require(Biostrings)  ## For reading FASTA via `readDNAStringSet'
require(Rsubread)    ## For Mapping and aligning
require(tools) ## For extracting basenames without ext

#' @importFrom rtracklayer import export
#' @importFrom Biostrings readDNAStringSet writeXStringSet width
#' @importFrom Rsubread buildindex align featureCounts
#' @importFrom tools file_path_sans_ext
#' @importFrom utils download.file read.table write.table
NULL
#> NULL

#' @title Get the URLs for the genome files
#' @description Dynamically generates URLs for the FASTA and GTF files of a
#'   given genome.
#' @param species String. Currently either `mus_musculus', `homo_sapiens', or
#'   `danio_rerio'. Default is `mus_musculus'.
#' @param build String. The build for the species of interest. If \code{NULL},
#'   it will use the latest build which is `GRCm39' for `mus_musculus' and
#'   `GRCh38' for `homo_sapiens'.
#' @param release String. The release version for that build. If \code{NULL}, it
#'   will use the latest release "105".
#' @param ensembl_base_url String depicting the base url of the Ensembl
#'   releases. Default is "http://ftp.ensembl.org/pub".
#' @return A list of two components. The first component \code{gtf} contains the
#'   URL of the GTF file on server. The second component \code{fasta} contains
#'   the URL of the FASTA primary assembly genome on server.
get_genome_urls <- function(species = "mus_musculus",
                            build = NULL, release = NULL,
                            ensembl_base_url = "http://ftp.ensembl.org/pub") {

  ## If other species are desired, they should be overridden directly in the
  ## urls_override parameter in the `get_genome_files' function.
  ensembl_species_map <- list(
    mus_musculus = list(build = "GRCm39", release = "105"),
    homo_sapiens = list(build = "GRCh38", release = "105"),
    danio_rerio = list(build = "GRCz11", release = "105")
  )

  generate_gtf_url <- function(species, build, release) {
    upper_species <- paste0(
      toupper(substring(species, 1, 1)),
      substring(species, 2, nchar(species))
    )
    paste0(
      ensembl_base_url, "/release-", release, "/gtf/",
      species, "/", upper_species, ".", build, ".", release, ".gtf.gz"
    )
  }

  generate_fasta_url <- function(species, build, release) {
    upper_species <- paste0(
      toupper(substring(species, 1, 1)),
      substring(species, 2, nchar(species))
    )
    paste0(
      ensembl_base_url, "/release-", release, "/fasta/",
      species, "/dna/", upper_species, ".", build,
      ".dna.primary_assembly.fa.gz"
    )
  }

  if (is.null(build)) {
    build <- ensembl_species_map[[species]]$build
  }
  if (is.null(release)) {
    release <- ensembl_species_map[[species]]$release
  }

  gtf_url <- generate_gtf_url(species, build, release)
  fasta_url <- generate_fasta_url(species, build, release)
  return(list(gtf = gtf_url, fasta = fasta_url))
}


#' @title Compute Unix checksums for downloaded genome files
#'
#' @description This function checks for existing ".sum" files which are
#'   prefixed to downloaded files and if not present, downloads them. It then
#'   computes the checksum locally and checks it with the downloaded.
#' @param url URL to remote genome file (.fa.gz or .gtf.gz)
#' @param outfile Local file where there genome file exists.
#' @return logical. If the file is consistent with that on the remote server.
check_sum_matches <- function(url, outfile) {
    parse_check_sum <- function(str) {
        res_split <- unlist(strsplit(str, split = " "))
        return(list(
            sum1 = as.integer(res_split[[1]]),
            sum2 = as.integer(res_split[[2]]),
            name = res_split[[3]]
        ))
    }

    get_check_sum <- function(url, outfile) {
        checksums_gtf <- paste0(dirname(url), "/CHECKSUMS")
        sumfile <- paste0(outfile, ".sum") ## this is not md5, it is unix "sum"
        if (file.exists(sumfile)) {
            message("Checksum: ", sumfile, " already exists")
            res <- readLines(sumfile)
        } else {
            message("Downloading: ", sumfile)
            download.file(checksums_gtf, sumfile)
            tab <- readLines(sumfile)
            res <- grep(basename(outfile), tab, value = TRUE)
            if (length(res) > 1) {
                message(res)
                stop("Multiple matches found")
            }
            message("Writing:", sumfile)
            writeLines(res, con = sumfile)
        }
        return(parse_check_sum(res))
    }

    remote_sum <- get_check_sum(url, outfile)
    local_sum <- parse_check_sum(
        system(paste0("sum ", outfile), intern = TRUE)
    )
    all(
        remote_sum$sum1 == local_sum$sum1,
        remote_sum$sum2 == local_sum$sum2
    )
}

#' @title Get the Genome files
#' @description Retrieve the local genome FASTA and GTF locations, and if they
#'   do not yet exist, download them.
#' @param species String. Currently either `mus_musculus', `homo_sapiens', or
#'   `danio_rerio'. Default is `mus_musculus'.
#' @param output_folder String of characters representing the directory to store
#'   the genome files. Will be created if non existent
#' @param download_timeout Positive integer number. If the downloads timeout,
#'   increase this to 10000. Default is 1000
#' @param check_sums logical. Perform consistency checks on downloaded
#'   files. Default is \code{TRUE}
#' @param build String. The build for the species of interest. If \code{NULL},
#'   it will use the latest build which is `GRCm39' for `mus_musculus' and
#'   `GRCh38' for `homo_sapiens'.
#' @param release String. The release version for that build. If \code{NULL}, it
#'   will use the latest release "105".
#' @param urls_override A list of strings with two components: `gtf' an URL for
#'   the Ensembl GTF file, `fasta' an URL for the Ensemble FASTA file. This
#'   parameter overrides the species, build, and release parameters. Default
#'   value is \code{NULL}.
#' @return A list of two components. The first component \code{gtf} contains the
#'   location of the GTF file on disk. The second component \code{fasta}
#'   contains the location of the FASTA primary assembly genome on disk.
#' @examples
#' \donttest{
#' mus_musc = get_genome_files(
#'     species = "mus_musculus",
#'     output_folder = "/genomes/"
#' )
#' }
#' @export
get_genome_files <- function(species = "mus_musculus",
                             output_folder = "genomes",
                             download_timeout = 1000,
                             check_sums = TRUE,
                             build = NULL,
                             release = NULL,
                             urls_override = NULL) {
    options(timeout = download_timeout)
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }
    if (is.null(urls_override)) {
        urls <- get_genome_urls(species, build, release)
    } else {
        message("Using override URLs")
        urls <- urls_override
    }

    outfile <- list()
    outfile$gtf <- file.path(output_folder, basename(urls$gtf))
    outfile$fasta <- file.path(output_folder, basename(urls$fasta))

    if (file.exists(outfile$gtf)) {
        message("GTF file: ", outfile$gtf, " exists already")
    } else {
        message("Downloading: ", basename(outfile$gtf))
        download.file(urls$gtf, outfile$gtf, method = "wget")
    }
    if (!(check_sums && check_sum_matches(urls$gtf, outfile$gtf))) {
        stop("Checksum for", outfile$gtf, "failed")
    }

    if (file.exists(outfile$fasta)) {
        message("FASTA file: ", outfile$gtf, " exists already, skipping")
    } else {
        message("Downloading: ", basename(outfile$fasta))
        download.file(urls$fasta, outfile$fasta, method = "wget")
    }
    if (!(check_sums && check_sum_matches(urls$fasta, outfile$fasta))) {
        stop("Checksum for", outfile$fasta, "failed")
    }
    return(outfile)
}

#' @title Add new sequences to both GTF and FASTA files
#' @description Import GTF and FASTA files, and append the new custom sequences
#'   to both, exporting them back to file.
#' @param genome_files A list of two components. The first component \code{gtf}
#'   contains the location of the GTF file on disk. The second component
#'   \code{fasta} contains the location of the FASTA primary assembly genome on
#'   disk.
#' @param new_seqs_file String denoting the location of the the new FASTA
#'   sequences to be inserted on disk. Multiple sequences to be inserted should
#'   all be placed into this file.
#' @return A list of two components. The first component \code{gtf} contains the
#'   location of the new annotated GTF file on disk. The second component
#'   \code{fasta} contains the location of the new annotated FASTA primary
#'   assembly genome on disk.
#' @examples
#' \donttest{
#' mus_musc = get_genome_files(
#'     species = "mus_musculus",
#'     output_folder = "/genomes/"
#' )
#' new_genome_files = add_seqs_to_gtf_and_fasta(mus_musc, "~/mynewseqs.fasta")
#' }
#' @export
add_seqs_to_gtf_and_fasta <- function(genome_files, new_seqs_file) {
  message("Inserting Sequences")
  insert_seqs <- readDNAStringSet(new_seqs_file)
  ## Remove leading spaces in names
  names(insert_seqs) <- trimws(names(insert_seqs))

  ## Output files
  gtf_file_noext <- file_path_sans_ext(file_path_sans_ext(genome_files$gtf))
  new_gtf_file <- paste0(gtf_file_noext, "_with_",
                         paste0(names(insert_seqs), collapse = "-"), ".gtf.gz")

  fasta_file_noext <- file_path_sans_ext(file_path_sans_ext(genome_files$fasta))
  new_genome_seq_file <- paste0(fasta_file_noext, "_with_",
                                paste0(names(insert_seqs), collapse = "-"),
                                ".fa.gz")

  if (file.exists(new_gtf_file)) {
    message("GTF exists: '", new_gtf_file, "' doing nothing.")
  } else {
    ## GTF
    message("Importing GTF")
    tabb <- as.data.frame(import(genome_files$gtf))
    ## -- Filter with CellRanger attributes
    ## Use these biotypes, as suggested by CellRanger
    keep_features <- tabb$gene_biotype %in%
      c(
        "protein_coding", "lincRNA", "antisense", "IG_LV_gene",
        "IG_V_gene", "IG_V_pseudogene", "IG_D_gene", "IG_J_gene",
        "IG_J_pseudogene", "IG_C_gene", "IG_C_pseudogene",
        "TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "TR_J_gene",
        "TR_J_pseudogene", "TR_C_gene"
      )
    message("GTF features before:                ", nrow(tabb))
    tabb <- tabb[keep_features, ]
    message("         - after biotype filter:  ", nrow(tabb))
    keep_genes <- !(is.na(tabb$gene_name) | is.null(tabb$gene_name))
    tabb <- tabb[keep_genes, ]
    message("         - after genename filter: ", nrow(tabb))

    message("Appending Sequences to GTF:")
    new_gtf_entry <- function(dnaseq_entry) {
      nam <- names(dnaseq_entry)
      len <- width(dnaseq_entry)
      src <- "ensembl"
      biotype <- "protein_coding"
      ver <- 1

      return(c(
        seqnames = nam, start = 1, end = len, width = len, strand = "+",
        source = src, type = "exon", score = NA, phase = NA,
        gene_id = nam, gene_version = ver, gene_name = nam,
        gene_source = src, gene_biotype = biotype, transcript_id = nam,
        transcript_version = ver, transcript_name = nam,
        transcript_source = src, transcript_biotype = biotype, tag = NA,
        transcript_support_level = NA, exon_number = 1, exon_id = nam,
        exon_version = ver, ccds_id = NA, protein_id = NA,
        protein_version = NA
      ))
    }

    add_lines_to_gtf <- function(gtf_table, insert_seqs) {
      levels(gtf_table$seqnames) <- c(
        levels(gtf_table$seqnames),
        names(insert_seqs)
      )
      for (ind in seq_along(names(insert_seqs))) {
        message("         - ", names(insert_seqs)[ind])
        gtf_table <- rbind(gtf_table, new_gtf_entry(insert_seqs[ind]))
      }
      return(gtf_table)
    }

    tab2 <- add_lines_to_gtf(tabb, insert_seqs)

    message("Saving new GTF file to: ", new_gtf_file)
    gtf_gz <- gzfile(new_gtf_file, "w")
    export(tab2, con = gtf_gz)
    close(gtf_gz)
  }

  if (file.exists(new_genome_seq_file)) {
    message("FASTA exists: '", new_genome_seq_file, "' doing nothing.")
  } else {
    ## FASTA
    message("Reading Genome FASTA")
    gen_seq <- readDNAStringSet(genome_files$fasta)
    message("Appending Sequences to FASTA")
    new_genome_seq <- c(gen_seq, insert_seqs)
    message("Saving new FASTA file to ", new_genome_seq_file)
    writeXStringSet(new_genome_seq, new_genome_seq_file,
                    append = FALSE, compress = TRUE)
  }
  return(list(gtf = new_gtf_file, fasta = new_genome_seq_file))
}

#' @title Prime FASTQ files for downstream
#' @description Discover FASTQ files either as single-end or paired-end.  Test
#'   them for consistency if paired-end, and prepare their BAM output names.
#' @param indir String denoting where the input directory of FASTQ files are.
#' @param r1_ending Suffix string that depicts the ending of the R1 reads.
#' @param r2_ending Suffix string the depicts the ending of the R2 reads. If
#'   \code{NULL}, then the FASTQ files are assumed to be single-end, and only
#'   the \code{r1_ending} reads are used for discovery. Default value is
#'   \code{NULL}.
#' @param align_end Suffix string that will be used as a suffix to name the BAM
#'   files.
#' @return A list of three components: reads1 contains a vector of R1 read
#'   filepaths, reads2 contains a vector of R2 read filepaths, and align_base
#'   contains the basename of the BAM files (no path).
#' @examples
#' \donttest{
#' dir_lists = list(
#'     index = index.dir,
#'     fastq = "1_fastqs",
#'     align = "2_bams",
#'     count = "3_counts"
#' )
#' read_lists = prime_fastq_files(
#'     dir_lists$fastq, "_1.fq.gz", "_2.fq.gz",
#'     "-align.bam"
#' )
#' }
#' @export
prime_fastq_files <- function(indir, r1_ending, r2_ending = NULL,
                              align_end = "align.bam") {
    reads_r1 <- list.files(path = indir, pattern = paste0("*", r1_ending))
    if (is.null(r2_ending)) {
        r2_ending <- ""
        reads_r2 <- NULL
    } else {
        reads_r2 <- list.files(path = indir, pattern = paste0("*", r2_ending))
    }
    readfile1 <- file.path(indir, reads_r1)
    readfile2 <- {
        if (is.null(reads_r2)) {
            message("We assume that these are single-end reads")
            c()
        } else {
            file.path(indir, reads_r2)
        }
    }
    if (!(all(file.exists(c(readfile1, readfile2))))) {
        stop("Not all files exist...")
    }
    if (!(all(sub(r1_ending, "", readfile1) == sub(r2_ending, "", readfile2)))) {
        stop("Could not match all R1 files to R2 files")
    }

    align_files_base <- basename(paste0(sub(r1_ending, "", readfile1), align_end))

    return(list(
        reads1 = readfile1, reads2 = readfile2,
        align_base = align_files_base
    ))
}


#' @title Generate or Retrieve Subread Index
#' @description Build a Subread index at the directory location of your genome
#'   FASTA file, or if it already exists, retrieve it.
#' @param genome_fasta String depicting the filepath of the genome FASTA file.
#' @return String depicting the filepath of the Subread index.
#' @examples
#' \donttest{
#' source("custom_genome.R")
#' new_genome_files = add_seqs_to_gtf_and_fasta(homo_sap, "GFP.fasta")
#' index_dir = retrieve_index(new_genome_files$fasta)
#' dir_lists = list(
#'     index = index_dir, fastq = "1_fastqs", align = "2_bams",
#'     count = "3_counts"
#' )
#' }
#' @export
retrieve_index <- function(genome_fasta) {
  index_dir <- paste0(file_path_sans_ext(file_path_sans_ext(genome_fasta)),
                      "-subread-index")
  if (dir.exists(index_dir)) {
    message("Index already built: ", index_dir)
  } else {
    dir.create(index_dir)
    message("Building Index at: ", index_dir)
    buildindex(basename = file.path(index_dir, "reference_index"),
               reference = genome_fasta)
  }
  return(index_dir)
}



#' @title Perform alignment of FASTQ files against Reference via Subread
#' @description Align FASTQ reads into BAM files via Subreads' `align' function.
#' @param dir_lists A list of directories, with two mandatory components:
#'   `index', containing the location of the subread index; and `align', the
#'   output directory of the BAM alignments.
#' @param read_lists A list of vectors, with three mandatory components:
#'   `reads1`, unnamed vector of FASTQ R1 reads; `reads2', unnamed vector of
#'   FASTQ R2 reads; `align_base', unnamed vector of future BAM file locations.
#' @param type String depicting input FASTQ type. Either "rna" or "dna". Default
#'   is "rna".
#' @param nthreads Positive integer for the number of threads to use. Default is
#'   30.
#' @return Void function. Output files are created at the location of the
#'   `dir_lists$align'.
#' @examples
#' \donttest{
#' read_lists = prime_fastq_files(dir_lists$fastq, "R1.fastq.gz", "R2.fastq.gz")
#' perform_alignment(dir_lists, read_lists)
#' }
#' @export
perform_alignment <- function(dir_lists, read_lists,
                              type = "rna", nthreads = 30) {

  stopifnot(c("index", "align") %in% names(dir_lists))
  stopifnot(c("reads1", "reads2", "align_base") %in% names(dir_lists))

  dir.create(dir_lists$index, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_lists$align, recursive = TRUE, showWarnings = FALSE)

  align(index = file.path(dir_lists$index, "reference_index"),
        readfile1 = read_lists$reads1,
        readfile2 = read_lists$reads2,
        type = type, # or 0 for RNA, 1 for DNA or "dna"
        output_format = "BAM",
        output_file = file.path(dir_lists$align,
                                read_lists$align_base),
        phredOffset = 33,
        nthreads = nthreads)
}

#' @title Summarize Subread alignment statistics
#' @description Collate all individual BAM statistics into a single table, and
#'   write the summaries to file.
#' @param dir_align String depicting the alignment directory
#' @param read_align A vector of strings depicting the names of the alignment
#'   files.
#' @return A data frame depicting the total fragments and the mapped fragments,
#'   along with the percentage mapped.
#' @examples
#' \donttest{
#' dir_lists = list(
#'     index = index_dir,
#'     fastq = "1_fastqs_remote_mount",
#'     align = "2_bams",
#'     count = "3_counts"
#' )
#' perform_alignment(dir_lists, read_lists)
#' read_lists = prime_fastq_files(dir_lists$fastq, "_1.fq.gz", "_2.fq.gz")
#' summarize_alignment(dir_lists$align, read_lists$align)
#' }
#' @export
summarize_alignment <- function(dir_align, read_align) {
  bam_summary <- paste0(file.path(dir_align, read_align), ".summary")

  if (!(all(file.exists(bam_summary)))) {
    message(bam_summary)
    stop("Not all summary files could be found")
  }

  tab <- as.matrix(do.call(
    cbind,
    lapply(bam_summary, function(x) {
      y <- as.data.frame(t(read.table(x, row.names = 1, check.names = FALSE)))
      perc_mapped <- floor(100 * y$Mapped_fragments / y$Total_fragments)
      y$Percentage_Mapped <- perc_mapped
      return(t(y))
    })
  ))

  colnames(tab) <- sub(
    "[._+:]align$", "",
    file_path_sans_ext(read_align)
  )

  stat_summ_file <- file.path(dir_align, "stats_summary.tsv")
  write.table(tab, stat_summ_file, sep = "\t", quote = FALSE)
  message("Wrote: ", stat_summ_file)
  return(t(as.data.frame(tab[c(
    "Total_fragments", "Mapped_fragments",
    "Percentage_Mapped"
  ), ])))
}


#' @title Generate Count Matrix from Aligned BAM files
#' @description Perform feature counts on aligned data and produce a count
#'   matrix and matrix statistics.
#' @param dir_lists A list of two mandatory components: `align' containing the
#'   directory of the alignment files, `count' containing the directory of where
#'   the count matrix will be placed.
#' @param gtf_file String depicting the filepath of where the GTF annotation
#'   file is located.
#' @param attrtype String depicting the attribute type to use as rownames in the
#'   final count matrix. All keys in the 9th column of the GTF file are
#'   possible. Default is "gene_name".
#' @param featuretype String depicting the feature type to count reads from. All
#'   values in the 3rd column of the GTF file are possible, such as `CDS',
#'   `exon', `five_prime_utr', `three_prime_utr', `gene', `transcript',
#'   etc. Default is "exon".
#' @param isPairedEnd logical. Whether the reads are paired-end or single-end.
#' @param nthreads Positive integer. Number of threads to use for
#'   counting. Default is 30.
#' @return Void function. Output matrices and statistics are placed in the
#'   `dir_lists$count' directory.
#' @examples
#' \donttest{
#' dir_lists = list(index=index_dir, fastq="1_fastqs",
#'                  align="2_bams", count="3_counts")
#' read_lists = prime_fastq_files(dir_lists$fastq, "R1.fastq.gz", "R2.fastq.gz")
#' perform_alignment(dir_lists, read_lists)
#' do_feature_counts(dir_lists, new_genome_files$gtf,
#'     featuretype = "exon",
#'     isPairedEnd = !is.null(read_lists$reads2))
#' }
#' @export
generate_count_matrix <- function(dir_lists, gtf_file,
                                  attrtype = "gene_name", featuretype = "exon",
                                  isPairedEnd = TRUE, nthreads = 30) {
  stopifnot(c("count", "align") %in% names(dir_lists))
  dir.create(dir_lists$count, recursive = TRUE, showWarnings = FALSE)

  bamfiles <- paste0(dir_lists$align, "/",
                     list.files(dir_lists$align, pattern = "*align.bam$"))
  names(bamfiles) <- basename(file_path_sans_ext(file_path_sans_ext(bamfiles)))

  message("Counting ", featuretype, " threads=", nthreads)
  fc <- featureCounts(files = bamfiles,
                      nthreads = nthreads,
                      isPairedEnd = isPairedEnd,
                      annot.ext = gtf_file,
                      isGTFAnnotationFile = TRUE,
                      GTF.featureType = featuretype,
                      GTF.attrType = attrtype,
                      juncCounts = TRUE)

  count_table <- fc$counts
  count_stat <- fc$stat
  count_table <- count_table[order(rownames(count_table)), ]

  cleaner_sample_names <- function(cnames) {
    cln <- basename(file_path_sans_ext(file_path_sans_ext(cnames)))
    return(sub("[._+:]align$", "", cln))
  }

  colnames(count_table) <- cleaner_sample_names(colnames(count_table))
  colnames(count_stat) <- cleaner_sample_names(colnames(count_stat))

  out_stat <- paste0(file.path(dir_lists$count,  featuretype),
                     ".", attrtype, ".stats.tsv")

  out_count <- paste0(file.path(dir_lists$count, featuretype),
                      ".", attrtype, ".matrix.tsv")

  ## Write Stats
  write.table(count_stat, out_stat, sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)
  message("Stats written to: ", out_stat)

  ## Write Counts
  write.table(count_table, out_count, sep = "\t",
              quote = FALSE, col.names = NA, row.names = TRUE)
  message("Counts written to: ", out_count)
}
