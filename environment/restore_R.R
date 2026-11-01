#!/usr/bin/env Rscript
# Exact, project-local restore; run from the repository root with R 4.6.1.
# No project activation, .Rprofile edit, global library write or version upgrade.
args <- commandArgs(trailingOnly = TRUE)
unknown <- args[!grepl('^--library=', args) & !args %in% c('--source', '--offline')]
if (length(unknown)) stop('Unknown arguments: ', paste(unknown, collapse=', '))
if (as.character(getRversion()) != '4.6.1') {
  stop('This lock was verified with R 4.6.1. Install that R version separately; this script never changes R.')
}
root <- normalizePath(getwd(), mustWork=TRUE)
if (!file.exists('environment/R-packages.lock.tsv')) stop('Run from the repository root.')
libarg <- args[grepl('^--library=', args)]
if (length(libarg) > 1) stop('Specify --library at most once.')
lib <- if (length(libarg)) sub('^--library=', '', libarg) else 'environment/R-library'
# Resolve the parent first to prevent writing through a path outside the project.
if (!dir.exists(dirname(lib))) stop('The library parent directory must already exist.')
parent <- normalizePath(dirname(lib), mustWork=TRUE)
if (!startsWith(paste0(parent, '/'), paste0(root, '/'))) stop('Library must stay inside the project.')
if (dir.exists(lib)) {
  resolved <- normalizePath(lib, mustWork=TRUE)
  if (!startsWith(paste0(resolved, '/'), paste0(root, '/'))) stop('Library symlink points outside the project.')
}
dir.create(lib, recursive=FALSE, showWarnings=FALSE)
lib <- normalizePath(lib, mustWork=TRUE)
.libPaths(c(lib, .Library))
options(Ncpus=2L, timeout=120)
Sys.setenv(OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', VECLIB_MAXIMUM_THREADS='2', MAKEFLAGS='-j2')
d <- read.delim('environment/R-packages.lock.tsv', check.names=FALSE, stringsAsFactors=FALSE)
if (any(!grepl('^[0-9a-f]{64}$',d$SHA256))) stop('Missing or invalid SHA256 in lock.')
if (anyDuplicated(paste(d$Package,d$Type))) stop('Duplicate package/archive entries in lock.')
is_arm_mac <- identical(Sys.info()[['sysname']], 'Darwin') && grepl('aarch64|arm64', R.version$platform)
prefer_binary <- is_arm_mac && !('--source' %in% args)
install_order <- unique(d$Package)
for (package in install_order) {
  rows <- d[d$Package == package, , drop=FALSE]
  version <- unique(rows$Version)
  if (length(version) != 1) stop('Conflicting locked versions for ',package)
  existing <- installed.packages(lib.loc=lib)
  if (package %in% rownames(existing)) {
    if (existing[package,'Version'] != version) {
      stop('Existing ',package,' has a different version in ',lib,
           '. Use a new project-local library; no package was overwritten.')
    }
    cat('Already exact:',package,version,'\n')
    next
  }
  archive_type <- if (prefer_binary && 'mac.binary' %in% rows$Type) 'mac.binary' else 'source'
  row <- rows[rows$Type == archive_type, , drop=FALSE]
  if (nrow(row) != 1) stop('No exact archive for ',package,' type=',archive_type)
  archive <- row$CachePath
  if (!startsWith(archive,'environment/package-cache/') || grepl('\\.\\.',archive)) stop('Unsafe cache path.')
  if (!file.exists(archive)) {
    if ('--offline' %in% args) stop('Offline cache missing: ',archive)
    dir.create(dirname(archive),recursive=TRUE,showWarnings=FALSE)
    partial <- paste0(archive,'.partial')
    download.file(row$URL,partial,mode='wb',method='libcurl',quiet=FALSE)
    if (unname(tools::sha256sum(partial)) != row$SHA256) stop('Downloaded checksum mismatch: ',package)
    if (!file.rename(partial,archive)) stop('Could not finalize archive: ',archive)
  }
  if (unname(tools::sha256sum(archive)) != row$SHA256) stop('Cached checksum mismatch: ',archive)
  cat('Installing exact archive:',package,version,archive_type,'\n')
  install.packages(archive,repos=NULL,type=archive_type,lib=lib,dependencies=FALSE,Ncpus=2L)
  installed <- installed.packages(lib.loc=lib)
  if (!package %in% rownames(installed) || installed[package,'Version'] != version) stop('Exact installation failed: ',package)
}
# Verify all manager-recorded packages, including transitive dependencies.
lock <- jsonlite::fromJSON('environment/renv.lock', simplifyVector=FALSE)
for (package in names(lock$Packages)) {
  expected <- lock$Packages[[package]]$Version
  actual <- packageDescription(package,lib.loc=lib)$Version
  if (actual != expected) stop('renv/archive lock disagreement for ',package)
}
cat('PASS: exact R package restore; library=',lib,'\n',sep='')
cat('R=',as.character(getRversion()),'; platform=',R.version$platform,'\n',sep='')
