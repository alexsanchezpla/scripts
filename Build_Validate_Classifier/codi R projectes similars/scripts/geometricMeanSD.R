geomean <- function(x, na.rm = FALSE, trim = 0, ...)
{
  exp(mean(log(x, ...), na.rm = na.rm, trim = trim, ...))
}

geosd <- function(x, na.rm = FALSE, ...)
{
  exp(sd(log(x, ...), na.rm = na.rm, ...))
}