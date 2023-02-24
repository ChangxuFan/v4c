#!/opt/htcf/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/r-3.6.3-wjv45rejzld4tqrpekyen7mtmm2ko6qo/bin/Rscript
suppressMessages(library(dplyr))
suppressMessages(library(BiocParallel))
suppressMessages(library("argparse"))
options(stringsAsFactors = F)

parser <- argparse::ArgumentParser(add_help = T)
parser$add_argument("-g", "--genome", required = T)
parser$add_argument("-c", "--threads.sub", required = T, type = "integer")
parser$add_argument("-C", "--threads.master", required = F, default = 1, type = "integer")

parser$add_argument("-d", "--dir", required = F, default="./", nargs = "+")

parser$add_argument("-s", "--sam", required = F, default = NULL, nargs= "+")
parser$add_argument("--subsample", required = F, default = NULL, type = "double")
parser$add_argument("--no_overwrite", required = F, action = "store_true", default = FALSE)

parser$add_argument("-t", "--txt", required = F, default = NULL, nargs = "+")
parser$add_argument("--htcf", required = F, action = "store_true", default = FALSE)
parser$add_argument("--slurm", required = F, action = "store_true", default = FALSE)

args <- parser$parse_args()


if (args$slurm == F) {
    print("using local parallel")
    par.param <- MulticoreParam(workers = args$threads.master, stop.on.error = F)
    if (args$htcf == F)
      samtools <- "/bar/cfan/anaconda2/envs/jupyter/bin/samtools"
    else
      samtools <- "/opt/htcf/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/samtools-1.12-awrseggdtfzhpv5y5evs43lz3b4eq472/bin/samtools"
  } else {
    print("using slurm")
    # stop("the slurm part has been working poorly")
    samtools <- "/opt/htcf/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/samtools-1.12-awrseggdtfzhpv5y5evs43lz3b4eq472/bin/samtools"
    TMPL.SLURM <- "/home/fanc/scripts/batchtools/templates/R3.6.3.tmpl"
    slurm.walltime = 172800
    slurm.memory = 4096

    slurm.dir <- paste0("./slurm",
     Sys.time() %>% as.character() %>% sub(" ", ":", .))
    # system(paste0("rm -rf ", slurm.dir))
    log.dir <- paste0(slurm.dir, "/log.dir/")
    result.dir <- paste0(slurm.dir, "/result.dir/")
    system(paste0("mkdir -p ", log.dir, " ", result.dir))
    
    par.param <- BatchtoolsParam(workers = args$threads.master, cluster = "slurm",
                                 registryargs = batchtoolsRegistryargs(
                                   file.dir = paste0(slurm.dir, "/file.dir/"),
                                   work.dir = "./",
                                   # packages = "BiocParallel",
                                   source = "~/R_packages/common/base.R"
                                 ),
                                 resources = list(walltime = slurm.walltime, memory = slurm.memory,
                                                  ncpus = args$threads.sub + 1),
                                 template = TMPL.SLURM,
                                 log = T, logdir = log.dir,
                                 resultdir = result.dir,
                                 stop.on.error = F)
  }


if (is.null(args$sam)) {
  juicer.sam <- Sys.glob(paste0(args$dir, "/splits/*.gz.sam"))
  if (!is.null(args$subsample)) {
    juicer.sam <- bptry({
      bplapply(juicer.sam, function(x) {
              juicer.sam <- utilsFanc::insert.name.before.ext(x, insert = args$subsample, "_")
              if (! file.exists(juicer.sam) || args$no_overwrite == F) {
                cmd <- paste0(samtools, " view -h -o ", juicer.sam, 
                  " -@ ", args$threads.sub,
                  " -s ", 42 + args$subsample, " ", x)
                print(cmd)
                system(cmd)
              }
              return(juicer.sam)
            }, BPPARAM = par.param) %>% unlist()
      })

  } 
	
}

# merge juicer.sam files located in the same directory
juicer.sam <- juicer.sam %>% split(., dirname(.)) %>%
  bplapply(function(x) {
      out.sam <- paste0(dirname(x[1]) , "/", 
        basename(x[1] %>% sub(".+\\.fastq", "merged.fastq", .)))
      cmd <- paste0(samtools, " merge -f -@ ", args$threads.sub, 
        " ", out.sam, " ", paste0(x, collapse = " "))
      print(cmd); system(cmd)
      return(out.sam)
    }, BPPARAM = par.param) %>% `names<-`(NULL) %>% unlist()

if (is.null(args$txt))
	nodups.txt <- Sys.glob(paste0(args$dir, "/aligned/merged_nodups.txt"))

print(juicer.sam)
print(nodups.txt)
# if (length(juicer.sam) != 1) {
# 	stop("length(juicer.sam) != 1")
# }

# if (length(nodups.txt) != 1) {
# 	stop("length(nodups.txt) != 1")
# }

if (length(juicer.sam) != length(nodups.txt)) {
	stop("length(juicer.sam) != length(nodups.txt)")
}



bptry({
	bplapply(1:length(juicer.sam), function(i) {
			v4c::juicer.txt.filter.ez(rname.file = NULL, juicer.sam = juicer.sam[i], juicer.txt = nodups.txt[i], 
                          mapq = 0, genome = args$genome, threads = args$threads.sub, run = T, samtools = samtools,
                          reverse.grep = T)
   #Sys.sleep(30)
		}, BPPARAM = par.param)
	})

