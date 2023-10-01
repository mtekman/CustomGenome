####################################################################
# required packages
####################################################################
require(rtracklayer) ## For reading GTF via `import'
require(Biostrings)  ## For reading FASTA via `DNASeqFromSet'
require(Rsubread)    ## For Mapping and aligning
require(tools) ## For extracting basenames without ext


#' @title Get the URLs for the genome files
#' @description Dynamically generates URLs for the FASTA and GTF files of a
#'   given genome.
#' @param species String. Currently either `mus_musculus' or
#'   `homo_sapiens'. Default is `mus_musculus'.
#' @param build String. The build for the species of interest. If \code{NULL},
#'   it will use the latest build which is `GRCm39' for `mus_musculus' and
#'   `GRCh38' for `homo_sapiens'.
#' @param release String. The release version for that build. If \code{NULL}, it
#'   will use the latest release "105".
#' @return A list of two components. The first component \code{gtf} contains the
#'   URL of the GTF file on server. The second component \code{fasta} contains
#'   the URL of the FASTA primary assembly genome on server.
get_genome_urls <- function(species = "mus_musculus", build = NULL,
                            release = NULL) {

  ensembl_base_url <- "http://ftp.ensembl.org/pub"
  species_map <- list(
    mus_musculus = list(build = "GRCm39", release = "105"),
    homo_sapiens = list(build = "GRCh38", release = "105")
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
      build <- species_map[[species]]$build
  }
  if (is.null(release)) {
      release <- species_map[[species]]$release
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
#' @param species String. Currently either `mus_musculus' or
#'   `homo_sapiens'. Default is `mus_musculus'.
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
#' @return A list of two components. The first component \code{gtf} contains the
#'   location of the GTF file on disk. The second component \code{fasta}
#'   contains the location of the FASTA primary assembly genome on disk.
#' @examples
#' mus_musc = getGenomeFiles(species="mus_musculus", output_folder="/genomes/")
#' @export
get_genome_files <- function(species = "mus_musculus",
                             output_folder = "genomes",
                             download_timeout = 1000,
                             check_sums = TRUE,
                             build = NULL,
                             release = NULL) {
  options(timeout = download_timeout)
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  urls <- get_genome_urls(species, build, release)

  outfile <- list()
  outfile$gtf <- file.path(output_folder, basename(urls$gtf))
  outfile$fasta <- file.path(output_folder, basename(urls$fasta))

  if (file.exists(outfile$gtf)) {
    message("GTF file: ", outfile$gtf, " exists already")
  } else {
    message("Downloading: ", basename(outfile$gtf))
    download.file(urls$gtf, outfile$gtf, method="wget")
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
#' mus_musc = getGenomeFiles(species="mus_musculus", output_folder="/genomes/")
#' new_genome_files = add_seqs_to_gtf_and_fasta(mus_musc, "~/mynewseqs.fasta")
#' @export
add_seqs_to_gtf_and_fasta <- function(genome_files, new_seqs_file) {

  message("Insert Seqs")
  insert_seqs <- readDNAStringSet(new_seqs_file)
  ## Remove leading spaces in names
  names(insert_seqs) <- trimws(names(insert_seqs))

  ## Output files
  gtf_file_noext <- file_path_sans_ext(file_path_sans_ext(genome_files$gtf))
  new_gtf_file <- paste0(gtf_file_noext, "_with_",
                         paste0(names(insert_seqs),
                                collapse = "-"), ".gtf.gz")

  fasta_file_noext <- file_path_sans_ext(file_path_sans_ext(genome_files$fasta))
  new_genome_seq_file <- paste0(fasta_file_noext, "_with_",
                                paste0(names(insert_seqs),
                                       collapse = "-"), ".fa.gz")

  if (file.exists(new_gtf_file)) {
    message("GTF exists: '", new_gtf_file, "' doing nothing.")
  } else {
    ## GTF
    message("Importing GTF")
    tabb <- as.data.frame(import(genome_files$gtf))
    ## -- Filter with CellRanger attributes
    ## Use these biotypes, as suggested by CellRanger
    keep_features <- tabb$gene_biotype %in%
      c("protein_coding", "lincRNA", "antisense", "IG_LV_gene",
        "IG_V_gene", "IG_V_pseudogene", "IG_D_gene", "IG_J_gene",
        "IG_J_pseudogene", "IG_C_gene", "IG_C_pseudogene",
        "TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "TR_J_gene",
        "TR_J_pseudogene", "TR_C_gene")
    message("GTF filter before:                ", nrow(tabb))
    tabb <- tabb[keep_features, ]      
    message("GTF filter after biotype filter:  ", nrow(tabb))
    keep_genes <- !(is.na(tabb$gene_name) | 
                    is.null(tabb$gene_name))
    tabb <- tabb[keep_genes, ]
    message("GTF filter after genename filter: ", nrow(tabb))
      
    message("Appending Sequences to GTF")
    new_gtf_entry <- function(dnaseq_entry) {
      nam <- names(dnaseq_entry)
      len <- width(dnaseq_entry)
      src <- "ensembl"
      biotype <- "protein_coding"
      ver <- 1

      return(c(seqnames = nam, start = 1, end = len, width = len, strand = "+",
               source = src, type = "exon", score = NA, phase = NA,
               gene_id = nam, gene_version = ver, gene_name = nam,
               gene_source = src, gene_biotype = biotype, transcript_id = nam,
               transcript_version = ver, transcript_name = nam,
               transcript_source = src, transcript_biotype = biotype, tag = NA,
               transcript_support_level = NA, exon_number = 1, exon_id = nam,
               exon_version = ver, ccds_id = NA, protein_id = NA,
               protein_version = NA))
    }

    add_lines_to_gtf <- function(gtf_table, insert_seqs) {
      levels(gtf_table$seqnames) <- c(levels(gtf_table$seqnames),
                                      names(insert_seqs))
      for (ind in seq_along(names(insert_seqs))) {
        message(names(insert_seqs)[ind])
        gtf_table <- rbind(gtf_table, new_gtf_entry(insert_seqs[ind]))
      }
      return(gtf_table)
    }

    tab2 <- add_lines_to_gtf(tabb, insert_seqs)
      
    message("Saving new GTF file to ", new_gtf_file)
    gtf_gz = gzfile(new_gtf_file, "w")
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
    writeXStringSet(new_genome_seq,
                    new_genome_seq_file, append = FALSE, compress = TRUE)
  }

  return(list(gtf = new_gtf_file, fasta = new_genome_seq_file))
}


prime_fastq_files <- function(indir, r1_ending, r2_ending, 
                              align_end="align.bam") {
    reads_R1 <- list.files(path = indir, pattern = paste0("*", r1_ending))
    if (is.null(r2_ending)){
        r2_ending <- ""
        reads_R2 <- NULL
    } else {
        reads_R2 <- list.files(path = indir, pattern = paste0("*", r2_ending))
    }    
    readfile1 = file.path(indir, reads_R1)
    readfile2 = {
        if (is.null(reads_R2)){
            message("We assume that these are single-end reads")
            c()
        } else {
            file.path(indir, reads_R2)
        }
    }
    if (!(all(file.exists(c(readfile1, readfile2))))) {
        stop("Not all files exist...")
    }
    if (!(all(sub(r1_ending, "", readfile1) == 
              sub(r2_ending, "", readfile2)))) {
        stop("Could not match all R1 files to R2 files") 
    }

    align_files_base <- basename(
        paste0(sub(r1_ending, "", readfile1), align_end))
    
    return(list(reads1 = readfile1, reads2 = readfile2, 
                align_base = align_files_base))
}

perform_alignment <- function(dir_lists, read_lists, nthreads=30, type="rna"){
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

summarize_alignment <- function(dir_align, read_lists){
    bam_summary = paste0(
        file.path(dir_lists$align,
                  read_lists$align_base), ".summary")

    if (!(all(file.exists(bam_summary)))){
        message(bam_summary)
        stop("Not all summary files could be found")
    }

    tab <- as.matrix(do.call(
            cbind, 
            lapply(bam_summary, function(x){
                y = as.data.frame(t(
                    read.table(x,row.names=1,check.names=FALSE)
                ))
                perc_mapped = floor(100 * y$Mapped_fragments / y$Total_fragments)
                y$Percentage_Mapped = perc_mapped
                return(t(y))
            })
    ))
    colnames(tab) <- sub("[._+:]align$", "",
                        file_path_sans_ext(read_lists$align_base))

    stat_summ_file <- file.path(dir_align, "stats_summary.tsv")
    write.table(tab,
                stat_summ_file, 
                sep="\t", quote=F)
    message("Wrote: ", stat_summ_file)
    return(t(as.data.frame(tab[c("Total_fragments", 
                               "Mapped_fragments", "Percentage_Mapped"),])))  
}


do_feature_counts <- function(dir_lists, gtf_file,
                              attrtype = "gene_name", 
                              featuretype = "exon", 
                              isPairedEnd=TRUE, nthreads=30){
    dir.create(dir_lists$count, recursive = TRUE, showWarnings = FALSE)
    
    bamfiles <- paste0(dir_lists$align, "/",
                       list.files(dir_lists$align, pattern="*align.bam$"))
    names(bamfiles) <- basename(file_path_sans_ext(file_path_sans_ext(bamfiles)))

    message("Counting ", featuretype, " threads=", nthreads)
    fc <- featureCounts(files = bamfiles,
                        nthreads = nthreads,
                        isPairedEnd = isPairedEnd,
                        annot.ext = gtf_file,
                        isGTFAnnotationFile = TRUE,
                        GTF.featureType = featuretype,
                        GTF.attrType = attrtype,
                        juncCounts = T)

    count_table <- fc$counts
    count_stat <- fc$stat
    count_table <- count_table[order(rownames(count_table)), ] 

    cleaner_sample_names <- function(cnames){
        cln = basename(file_path_sans_ext(file_path_sans_ext(cnames)))
        return(sub("[._+:]align$", "", cln))
    }
    
    colnames(count_table) <- cleaner_sample_names(colnames(count_table))
    colnames(count_stat) <- cleaner_sample_names(colnames(count_stat))

    out_stat <- paste0(file.path(dir_lists$count,  featuretype), 
                       ".", attrtype, ".stats.tsv")

    out_count <- paste0(file.path(dir_lists$count, featuretype), 
                       ".", attrtype, ".matrix.tsv")
    
    ## Write Stats
    write.table(count_stat, out_stat,
        sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    message("Stats written to: ", out_stat)

    ## Write Counts for each BAM
    ##zzz <- lapply(colnames(count_table), function(nme){    
    ##    out_table <- paste0(
    ##        file.path(out_table_dir, nme), attrtype, ".tab")
    ##    write.table(count_table[,nme], out_table,
    ##                sep = "\t", quote = FALSE, 
    ##                col.names = FALSE)
    ##    message("Table written to: ", out_table)
    ##})

    ## Write Counts
    write.table(count_table, out_count,
                sep="\t", quote=FALSE, col.names=NA, row.names=TRUE)
    message("Counts written to: ", out_count)

    
}
