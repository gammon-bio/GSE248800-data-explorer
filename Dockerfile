# Dockerfile — self-contained image for the GSE248800 non-myofiber explorer.
#
# The app bundles ~675 MB of HDF5 expression data, so a container (data baked in,
# no upload limits) is the natural deployment. Runtime needs only the shiny stack
# + rhdf5 — fgsea/progeny/msigdbr are used only offline in data-raw/, never here.
#
# Build (from repo root):   docker build -t bryce-scrnaseq .
# Run:                      docker run --rm -p 3838:3838 bryce-scrnaseq
# Then open:                http://localhost:3838
#
# Deploys unchanged to Google Cloud Run, a VPS with Docker, or anywhere that runs
# a container. rocker/r-ver pins R 4.5.3 and installs Linux binaries from Posit
# Public Package Manager (P3M), so builds are fast and reproducible.

FROM rocker/r-ver:4.5.3

# System libraries needed to install/run the R package set. rhdf5 bundles its own
# HDF5 (Rhdf5lib), so no system libhdf5 is required.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev libpng-dev \
      libfontconfig1-dev libfreetype6-dev libjpeg-dev \
      libuv1 libicu-dev libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# CRAN runtime deps (P3M Linux binaries via rocker's default repo).
RUN R -q -e "install.packages(c( \
      'shiny','bslib','shinyWidgets','DT','plotly','ggplot2','dplyr','tidyr', \
      'viridis','scales','waiter','htmltools'))"

# rhdf5 (Bioconductor) — the on-disk expression store backend.
RUN R -q -e "install.packages('BiocManager'); BiocManager::install('rhdf5', update=FALSE, ask=FALSE)"

# Fail the build early if any runtime package is missing.
RUN R -q -e "pkgs<-c('shiny','bslib','shinyWidgets','DT','plotly','ggplot2','dplyr','tidyr','viridis','scales','waiter','htmltools','rhdf5'); m<-pkgs[!pkgs %in% rownames(installed.packages())]; if(length(m)) stop('missing: ', paste(m,collapse=', ')) else cat('all runtime packages present\n')"

# App + baked-in data (global.R uses app_data_dir='data' relative to the app dir).
# Copy the large data FIRST (changes rarely) as its own layer, then the small
# code files — so an app-code edit re-pushes only a few KB, not the ~700 MB
# expression layer.
WORKDIR /srv/app
COPY app/data/ /srv/app/data/
COPY app/R/ /srv/app/R/
COPY app/modules/ /srv/app/modules/
COPY app/www/ /srv/app/www/
COPY app/app.R app/global.R app/ui.R app/server.R /srv/app/

# Railway/Cloud Run inject $PORT and health-check against it; honour it, and
# fall back to 3838 for a plain `docker run`.
EXPOSE 3838
CMD ["R", "-q", "-e", "shiny::runApp('/srv/app', host='0.0.0.0', port=as.integer(Sys.getenv('PORT','3838')))"]
