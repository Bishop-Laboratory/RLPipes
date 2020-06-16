# RSeq
A biologist-friendly best-practices R-loop mapping pipeline 

## Quick Start
RSeq can be initialized through `R`, command line, or in a user-friendly web interface. 

### Installation
It is also highly recommended to install RSeq through `conda`:
``` 
conda install -c bioconda rseq
```
This will install the command line tool, web launcher, and `R` package at once.
It will also install `snakemake` and all other non-`R` dependencies.

### Using RSeq through a web interface
RSeq can be run as a web application using `R-Shiny.` To run this, type:
```
rseq launch
```
This will launch the web interface which is used for configuring and monitoring 
the RSeq pipeline. By default, the web interface will be instructed to listen for 
traffic on port `6644` coming from `localhost`. It can be accessed by navigating in a 
browser to `localhost:6644`. 

However, it is not very convenient to interact this way with RSeq if running it on a 
remote server as this would require opening a browser on the remote machine and forwarding
the graphic interface through `X11`. Instead, users can specify 

