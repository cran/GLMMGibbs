.First.lib <- function(lib, pkg)
{
   library.dynam("GLMMGibbs", pkg, lib)
   invisible()
}

