####################################################################
## required packages
####################################################################

#' @importFrom rtracklayer import export
#' @importFrom Biostrings readDNAStringSet writeXStringSet width
#' @importFrom GenomicRanges GRanges
#' @importFrom Rsubread buildindex align featureCounts
#' @importFrom tools file_path_sans_ext
#' @importFrom utils download.file read.table write.table capture.output
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd

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
#' @param fasta_type String. The type of genomic FASTA sequence to retrieve.
#'   Default is "dna.primary_assembly",
#' @param gtf_type String. The type of genomic GTF file to retrieve.
#'   Default is "gtf", Other examples are "chr.gtf" or "abinitio.gtf".
#' @param ensembl_base_url String depicting the base url of the Ensembl
#'   releases. Default is "http://ftp.ensembl.org/pub".
#' @return A list of two components. The first component \code{gtf} contains the
#'   URL of the GTF file on server. The second component \code{fasta} contains
#'   the URL of the FASTA primary assembly genome on server.
#' @examples
#' CustomGenome:::get_genome_urls(
#'     species = "mus_musculus",
#'     release = "104"
#' )
get_genome_urls <- function(species = "mus_musculus",
                            build = NULL, release = NULL,
                            fasta_type = "dna.primary_assembly",
                            gtf_type = "gtf",
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
            species, "/", upper_species, ".", build, ".",
            release, ".", gtf_type, ".gz"
        )
    }

    generate_fasta_url <- function(species, build, release) {
        upper_species <- paste0(
            toupper(substring(species, 1, 1)),
            substring(species, 2, nchar(species))
        )
        paste0(
            ensembl_base_url, "/release-", release, "/fasta/",
            species, "/dna/", upper_species, ".", build, ".",
            fasta_type, ".fa.gz"
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

#' @title Get the Genome files
#' @description Retrieve the local genome FASTA and GTF locations, and if they
#'   do not yet exist, download them.
#' @param species String. Currently either `mus_musculus', `homo_sapiens', or
#'   `danio_rerio'. Default is `mus_musculus'.
#' @param cache_folder String of characters representing the directory to store
#'   the genome files. Will be created if non existent
#' @param download_timeout Positive integer number. If the downloads timeout,
#'   increase this to 10000. Default is 1000
#' @param urls_override A list of strings with two components: `gtf' an URL for
#'   the Ensembl GTF file, `fasta' an URL for the Ensemble FASTA file. This
#'   parameter overrides the species, build, and release parameters. Default
#'   value is \code{NULL}.
#' @param ... Arguments to be passed onto \code{get_genome_urls}.
#' @return A list of two components. The first component \code{gtf} contains the
#'   location of the GTF file on disk. The second component \code{fasta}
#'   contains the location of the FASTA primary assembly genome on disk.
#' @examples
#' mus_musc = get_genome_files(
#'     species = "mus_musculus",
#'     cache_folder = tempdir(),
#'     gtf_type = "abinitio.gtf",
#'     download_timeout = 10000
#' )
#' @export
get_genome_files <- function(species = "mus_musculus",
                            cache_folder = "genomes",
                            download_timeout = 1000,
                            urls_override = NULL, ...) {
    options(timeout = download_timeout)
    if (!dir.exists(cache_folder)) {
        dir.create(cache_folder)
    }
    if (is.null(urls_override)) {
        urls <- get_genome_urls(...)
    } else {
        message("Using override URLs")
        urls <- urls_override
    }

    outfile <- list()
    base_gtf <- basename(urls$gtf)      ## Keys to access
    base_fasta <- basename(urls$fasta)  ## records

    bfc <- BiocFileCache(cache_folder, ask=FALSE)
    gtf_record <- bfcquery(bfc, base_gtf)
    fasta_record <- bfcquery(bfc, basename(urls$fasta))

    if (nrow(gtf_record) != 0) {
        message("GTF file: ", base_gtf, "\n          exists already")
        outfile$gtf <- gtf_record$rpath[1]
    } else {
        message("Downloading: ", base_gtf)
        outfile$gtf <- bfcadd(bfc, base_gtf, fpath=urls$gtf)
    }
    message("Using GTF file: ", outfile$gtf)

    if (nrow(fasta_record) != 0) {
        message("FASTA file: ", base_fasta, "\n            exists already")
        outfile$fasta <- fasta_record$rpath[1]
    } else {
        message("Downloading: ", base_fasta)
        outfile$fasta <- bfcadd(bfc, base_fasta, fpath=urls$fasta)
    }
    message("Using FASTA file: ", outfile$fasta)
    return(outfile)
}

#' @title Filter sequences in a GTF for quality control
#' @description Restrict available sequences to the same biotypes used
#'   by CellRanger
#'   \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr}
#' @param gtf_file String depicting the filename of a GTF file to modify
#' @param wanted_biotypes An unnamed vector of biotypes to keep.
#'    Default uses the full list of biotypes from CellRanger.
#' @return A GTFFile object with filtered sequences.
#' @examples
#' sub8 <- system.file("extdata", "subsample8.gtf.gz", package = "CustomGenome")
#' new_gtf <- CustomGenome:::qc_filter_lines_of_gtf(sub8)
qc_filter_lines_of_gtf <- function(gtf_file,
                                    wanted_biotypes = c(
                                        "protein_coding", "lincRNA", "antisense",
                                        "IG_LV_gene", "IG_V_gene",
                                        "IG_V_pseudogene", "IG_D_gene",
                                        "IG_J_gene", "IG_J_pseudogene",
                                        "IG_C_gene", "IG_C_pseudogene",
                                        "TR_V_gene", "TR_V_pseudogene",
                                        "TR_D_gene", "TR_J_gene",
                                        "TR_J_pseudogene", "TR_C_gene"
                                    )) {
    message("Importing GTF")
    tabb <- import(gtf_file) ## Do not coerce into data frame
    ## -- Filter with CellRanger attributes
    ## Use these biotypes, as suggested by CellRanger
    keep_features <- tabb$gene_biotype %in% wanted_biotypes

    message("GTF features before:                ", length(tabb))
    tabb <- tabb[keep_features, ]
    message("         - after biotype filter:  ", length(tabb))
    keep_genes <- !(is.na(tabb$gene_name) | is.null(tabb$gene_name))
    tabb <- tabb[keep_genes, ]
    message("         - after genename filter: ", length(tabb))
    return(tabb)
}

#' @title Append Sequences to GTF file
#' @description Insert all sequences as single Ensembl "protein_coding" exons
#' @param gtffile A GTFFile object (GRanges), explicitly not a data frame.
#' @param insert_seqs A DNAStringSet object of all sequences.
#' @param ftype String depicting the feature type to add. Can be gene,
#'    exon or any other feature type. Default: exon
#' @param src String depicting source of sequence. Default: ensembl.
#' @param biotype String depicting the biotype. Default: protein_coding.
#' @return A modified GTFFile object with appended sequences included.
#' @examples
#' sub8 <- system.file("extdata", "subsample8.gtf.gz", package = "CustomGenome")
#' user_sequences <- system.file("extdata", "user_sequences.fa",
#'                                 package = "CustomGenome")
#' library(rtracklayer)
#' library(Biostrings)
#' obj_gtf <- import(sub8)
#' obj_fas <- readDNAStringSet(user_sequences)
#' tail(CustomGenome:::add_lines_to_gtf(obj_gtf, obj_fas))
add_lines_to_gtf <- function(gtffile, insert_seqs,
                            ftype = "exon", src = "ensembl",
                            biotype = "protein_coding") {
    new_gtf_entry <- function(dnaseq_entry) {
        nam <- names(dnaseq_entry)
        len <- width(dnaseq_entry)
        src <- "ensembl"
        biotype <- "protein_coding"
        ver <- 1

        return(GRanges(
            paste0(nam, ":", 1, "-", len, ":+"),
            source = src, type = ftype,
            gene_id = nam, gene_version = ver, gene_name = nam,
            gene_source = src, gene_biotype = biotype, transcript_id = nam,
            transcript_version = ver, transcript_name = nam,
            transcript_source = src, transcript_biotype = biotype,
            exon_number = 1, exon_id = nam, exon_version = ver
        ))
    }
    message("Appending Sequences to GTF:")
    for (ind in seq_along(names(insert_seqs))) {
        message("    - ", names(insert_seqs[ind]))
        suppressWarnings(  ## A merge warning might confuse the user here.
            gtffile <- c(gtffile, new_gtf_entry(insert_seqs[ind]))
        )
    }
    return(gtffile)
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
#' myfasta <- system.file("extdata", "tiny.fa.gz", package = "CustomGenome")
#' sub8 <- system.file("extdata", "subsample8.gtf.gz", package = "CustomGenome")
#' user_sequences <- system.file("extdata", "user_sequences.fa",
#'                                 package="CustomGenome")
#' gfiles <- list(gtf = sub8, fasta = myfasta)
#' new_genome_files <- add_seqs_to_gtf_and_fasta(gfiles, user_sequences)
#' @export
add_seqs_to_gtf_and_fasta <- function(genome_files, new_seqs_file) {
    message("Inserting Sequences")
    insert_seqs <- readDNAStringSet(new_seqs_file)
    ## Remove leading spaces in names
    names(insert_seqs) <- trimws(names(insert_seqs))

    ## Output files
    gtf_file_noext <- file_path_sans_ext(file_path_sans_ext(genome_files$gtf))
    new_gtf_file <- paste0(gtf_file_noext, "_with_",
                            paste0(names(insert_seqs), collapse = "-"),
                            ".gtf.gz")

    fasta_file_noext <- file_path_sans_ext(file_path_sans_ext(
        genome_files$fasta))
    new_genome_seq_file <- paste0(fasta_file_noext, "_with_",
                                    paste0(names(insert_seqs), collapse = "-"),
                                    ".fa.gz")

    if (file.exists(new_gtf_file)) {
        message("GTF exists: '", new_gtf_file, "' doing nothing.")
    } else {
        ## GTF
        tabb <- qc_filter_lines_of_gtf(genome_files$gtf)
        tabb <- add_lines_to_gtf(tabb, insert_seqs)

        message("Saving new GTF file to: ", new_gtf_file)
        gtf_gz <- gzfile(new_gtf_file, "w")
        export(tabb, con = gtf_gz)
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
#' prime_fastq_files(
#'         system.file("extdata", package = "CustomGenome"),
#'         "e4.gtf.gz", "e8.gtf.gz", "e-align.bam"
#' )
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
    if (!(all(sub(r1_ending, "", readfile1) ==
                sub(r2_ending, "", readfile2)))) {
        stop("Could not match all R1 files to R2 files")
    }

    align_files_base <- basename(paste0(sub(r1_ending, "",
                                            readfile1), align_end))

    return(list(
        reads1 = readfile1, reads2 = readfile2,
        align_base = align_files_base
    ))
}

#' @title Generate or Retrieve Subread Index
#' @description Build a Subread index at the directory location of your genome
#'   FASTA file, or if it already exists, retrieve it.
#' @param genome_fasta String depicting the filepath of the genome FASTA file.
#' @param index_dir String depicting the filepath of where the index directory
#'   will be built. If \code{NULL}, then it will be in the same location as
#'   the genome_fasta with "-subread-index" appended.
#' @param ... Extra arguments to be passed into Subreads `buildindex' function.
#' @return String depicting the filepath of the Subread index.
#' @examples
#' tiny <- system.file("extdata", "tiny.fa.gz", package = "CustomGenome")
#' retrieve_index(tiny)
#' @export
retrieve_index <- function(genome_fasta, index_dir = NULL, ...) {
    if (is.null(index_dir)) {
        index_dir <- paste0(
            file_path_sans_ext(file_path_sans_ext(genome_fasta)),
            "-subread-index")
    }
    if (dir.exists(index_dir)) {
        message("Index already built: ", index_dir)
    } else {
        dir.create(index_dir)
        message("Building Index at: ", index_dir)
        buildindex(
            basename = file.path(index_dir, "reference_index"),
            reference = genome_fasta,
            ...
        )
    }
    return(file.path(index_dir, "reference_index"))
}


#' @title Validate the align Ellipsis Arguments for RSubread's align function
#' @description Convert ellipsis arguments into a list and validate their value.
#' @param ... An ellipsis argument that is passed in Rsubreads's `align'
#'   function. Typical values are `featuretype' (default: "exon"), `isPairedEnd'
#'   (default: TRUE), `nthreads' (default: 30).
#' @return A list containing validated and/or modified ellipsis arguments.
validate_align_arguments <- function(...) {
    ## Set default args for featurecounts
    ellips <- list(...)
    ellips$type <- ifelse(is.null(ellips$type),
                        "rna", ellips$type)
    ellips$phredOffset <- ifelse(is.null(ellips$phredOffset),
                                33, ellips$phredOffset)
    ellips$output_format <- ifelse(is.null(ellips$output_format),
                                    "BAM", ellips$output_format)
    return(ellips)
}

#' @title Perform alignment of FASTQ files against Reference via
#'   Subread
#' @description Align FASTQ reads into BAM files via Subreads' `align'
#'   function.
#' @param dir_lists A list of directories, with three mandatory
#'   components: `index', containing the location of the subread
#'   index; `align', the output directory of the BAM alignments; and
#'   `stats', directory for stats and logging files.
#' @param read_lists A list of vectors, with three mandatory
#'   components: `reads1`, unnamed vector of FASTQ R1 reads; `reads2',
#'   unnamed vector of FASTQ R2 reads; `align_base', unnamed vector of
#'   future BAM file locations.
#' @param nthreads Positive integer for the number of threads to use.
#'   Default is 8.
#' @param ... An ellipsis argument that is passed in Rsubreads's
#'   `align' function. Typical values are `type' (default: "rna"),
#'   `phredOffset' (default: 33), `output_format' (default: "BAM").
#' @return Void function. Output files are created at the location of
#'   the `dir_lists$align'.
#' @examples
#' test_dir <- tempdir()
#' tiny <- system.file("extdata", "tiny.fa.gz", package = "CustomGenome")
#' r1_fastq <- system.file("extdata", "r1.fq.gz", package = "CustomGenome")
#' r2_fastq <- system.file("extdata", "r2.fq.gz", package = "CustomGenome")
#' read_lists <- list(
#'     reads1 = c(r1_fastq), reads2 = c(r2_fastq),
#'     align_base = c("align.bam")
#' )
#' dir_lists <- list(
#'     index = retrieve_index(tiny, paste0(test_dir, "-subread-index")),
#'     align = paste0(test_dir, "-aligned"),
#'     count = paste0(test_dir, "-counts"),
#'     stats = paste0(test_dir, "-stats")
#' )
#' perform_alignment(dir_lists, read_lists, nthreads = 1)
#' @export
perform_alignment <- function(dir_lists, read_lists, nthreads = 8, ...) {
    stopifnot(c("index", "align", "stats") %in% names(dir_lists))
    stopifnot(c("reads1", "reads2", "align_base") %in% names(read_lists))

    dir.create(dir_lists$align, recursive = TRUE, showWarnings = FALSE)
    dir.create(dir_lists$stats, recursive = TRUE, showWarnings = FALSE)

    ellips <- validate_align_arguments(...)
    ## Add required arguments
    ellips$index <- dir_lists$index
    ellips$readfile1 <- read_lists$reads1
    ellips$readfile2 <- read_lists$reads2
    ellips$nthreads <- nthreads
    ellips$output_file <- file.path(
        dir_lists$align,
        read_lists$align_base
    )
    message("Alignment : threads=", ellips$nthreads)
    out_log <- file.path(dir_lists$stats, "align.log")
    std_out <- capture.output(
        do_align <- do.call(align, ellips),
        file = out_log,
        split = TRUE
    )
    message("Alignment : finished")
    return(do_align)
}

#' @title Summarize Subread alignment statistics
#' @description Collate all individual BAM statistics into a single
#'   table, and write the summaries to file.
#' @param dir_lists A list of directories, with two mandatory
#'   components: `align', containing the location of the BAM
#'   alignments; and `stats' indicating where to store the statistics.
#' @param read_align A vector of strings depicting the base names of
#'   the alignment files.
#' @return A data frame depicting the total fragments and the mapped
#'   fragments, along with the percentage mapped.
#' @export
summarize_alignment <- function(dir_lists, read_align) {
    stopifnot(c("align", "stats") %in% names(dir_lists))
    bam_summary <- paste0(file.path(dir_lists$align, read_align), ".summary")

    if (!(all(file.exists(bam_summary)))) {
        message(bam_summary)
        stop("Not all summary files could be found")
    }

    tab <- as.matrix(do.call(
        cbind,
        lapply(bam_summary, function(x) {
            y <- as.data.frame(t(read.table(x, row.names = 1,
                                            check.names = FALSE)))
            perc_mapped <- floor(100 * y$Mapped_fragments / y$Total_fragments)
            y$Percentage_Mapped <- perc_mapped
            return(t(y))
        })
    ))

    colnames(tab) <- sub(
        "[._+:]align$", "",
        file_path_sans_ext(read_align)
    )

    stat_summ_file <- file.path(dir_lists$stats, "align.stats.tsv")
    write.table(tab, stat_summ_file, sep = "\t", quote = FALSE)
    message("Wrote: ", stat_summ_file)
    tab_out <- as.data.frame(t(tab))[
        c("Total_fragments", "Mapped_fragments", "Percentage_Mapped")]
    return(tab_out)
}


#' @title Validate the featureCounts Ellipsis Arguments for RSubread's
#'   featureCounts function.
#' @description Convert ellipsis arguments into a list and validate
#'   their value.
#' @param ... An ellipsis argument that is passed in Rsubreads's
#'   `featureCounts' function. Typical values are `GTF.featuretype'
#'   (default: "exon"), `isPairedEnd' (default: TRUE).
#' @return A list containing validated and/or modified ellipsis
#'   arguments.
validate_fc_arguments <- function(...) {
    ## Set default args for featurecounts
    ellips <- list(...)
    ellips$GTF.attrType <- ifelse(is.null(ellips$GTF.attrType), "gene_name",
                                    ellips$GTF.attrType)
    ellips$GTF.featureType <- ifelse(is.null(ellips$GTF.featureType), "exon",
                                    ellips$GTF.featureType)
    ellips$isPairedEnd <- ifelse(is.null(ellips$isPairedEnd), TRUE,
                                ellips$isPairedEnd)
    ellips$juncCounts <- ifelse(is.null(ellips$juncCounts), TRUE,
                                ellips$juncCounts)
    return(ellips)
}

#' @title Generate Count Matrix from Aligned BAM files
#' @description Perform feature counts on aligned data and produce a
#'   count matrix and matrix statistics.
#' @param dir_lists A list of three mandatory components: `align'
#'   containing the directory of the alignment files, `count'
#'   containing the directory of where the count matrix will be
#'   placed, `stats' containing the directory of where the stats
#'   tables will be placed.
#' @param gtf_file String depicting the filepath of where the GTF
#'   annotation file is located.
#' @param bam_pattern Pattern string depicting a way to match the BAM
#'   files in the `align' directory. Default is "*.bam$".
#' @param nthreads Positive integer for the number of threads to use.
#'   Default is 8.
#' @param ... An ellipsis argument that is passed in Rsubreads's
#'   `featurecounts' function. Typical values are `featuretype'
#'   (default: "exon"), `isPairedEnd' (default: TRUE).
#' @return A list of two components: `count' for the location of the
#'   count matrix file, and `stat' for the location of the statistics
#'   file. Output matrices and statistics are placed in the
#'   `dir_lists$count' directory.
#' @examples
#' dir_lists <- list(align = system.file("extdata", package="CustomGenome"),
#'                     count = tempdir(), stats = tempdir())
#' gtf_file <- system.file("extdata", "tiny.gtf.gz", package="CustomGenome")
#' generate_count_matrix(dir_lists, gtf_file, bam_pattern="*test.bam$")
#' @export
generate_count_matrix <- function(dir_lists, gtf_file, bam_pattern = "*.bam$",
                                nthreads = 8, ...) {
    stopifnot(c("count", "align", "stats") %in% names(dir_lists))
    dir.create(dir_lists$count, recursive = TRUE, showWarnings = FALSE)
    found_bams <- list.files(dir_lists$align, pattern = bam_pattern)
    message("Bam files found with pattern `", bam_pattern, "':\n - ",
            paste0(found_bams, collapse = "\n - "))
    if (length(found_bams) == 0) {
        stop("No Bam files found. Perhaps adjust the `bam_pattern' arguments")
    }
    bamfiles <- paste0(dir_lists$align, "/", found_bams)
    names(bamfiles) <- basename(file_path_sans_ext(
        file_path_sans_ext(bamfiles)))
    ellips <- validate_fc_arguments(...)
    ellips$files <- bamfiles; ellips$annot.ext <- gtf_file
    ellips$isGTFAnnotationFile <- TRUE; ellips$nthreads <- nthreads
    message("Counting : ", ellips$GTF.featureType, " threads=", ellips$nthreads)
    out_log <- file.path(dir_lists$stats,
                        paste0("count.", ellips$GTF.featureType, ".",
                                ellips$GTF.attrType, ".log"))
    std_out <- capture.output(fc <- do.call(featureCounts, ellips),
                                file = out_log, split = TRUE)
    message("Counting : finished")
    count_table <- fc$counts;count_stat <- fc$stat;
    count_table <- count_table[order(rownames(count_table)), ]
    cleaner_sample_names <- function(cnames) {
        cln <- basename(file_path_sans_ext(file_path_sans_ext(cnames)))
        return(sub("[._+:]align$", "", cln))}
    ## In the case of a single BAM file, still give it a name
    if (is.null(colnames(count_table))) {
        count_table <- as.data.frame(count_table)
        colnames(count_table) <- names(bamfiles)}
    colnames(count_table) <- cleaner_sample_names(colnames(count_table))
    colnames(count_stat) <- cleaner_sample_names(colnames(count_stat))
    out_stat <- file.path(dir_lists$stats,
                            paste0("count.", ellips$GTF.featureType, ".",
                                    ellips$GTF.attrType, ".stats.tsv"))
    out_count <- file.path(dir_lists$count,
                            paste0(ellips$GTF.featureType, ".",
                                    ellips$GTF.attrType, ".matrix.tsv"))
    ## Write Stats
    write.table(count_stat, out_stat, sep = "\t",
                quote = FALSE, col.names = TRUE, row.names = FALSE)
    message("Stats written to: ", out_stat)
    ## Write Counts
    write.table(count_table, out_count, sep = "\t",
                quote = FALSE, col.names = NA, row.names = TRUE)
    message("Counts written to: ", out_count)
    return(list(count = out_count, stat = out_stat))
}
