.onLoad <- function(libname, pkgname)
{
  ns <- asNamespace(pkgname)
	file <- system.file("extdata", "Zmays.angsd.tilbkn.AGPv2.17.Rdata",
										 	package=pkgname, lib.loc=libname)
	objname <- sub(".Rdata$", "", basename(file))
	load(file, envir=ns)
	namespaceExport(ns, objname)
}
