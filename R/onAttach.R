.onAttach<-function(libname, pkgname){
  pkg.desc <- read.dcf(file = system.file("DESCRIPTION", package = pkgname, lib.loc = libname))
  rg.desc <- capture.output(write.dcf(pkg.desc))
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("--------------------------------------------------------\n", appendLF = FALSE)
  packageStartupMessage("> The ReGenesees package has been successfully loaded. <\n", appendLF = FALSE)
  packageStartupMessage("--------------------------------------------------------\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  packageStartupMessage("\n", appendLF = FALSE)
  for (line in rg.desc) {
         packageStartupMessage(line)
    }
  packageStartupMessage("\n", appendLF = FALSE)
}
