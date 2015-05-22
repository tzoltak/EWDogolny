#' @title Predykcje z modeli mieszanych efektow
#' @description
#' Funkcja zwraca predykcje realizacji ef. los. i związane z nimi warunkowe (a'posteriori)
#' odchylenia standardowe dla poszczególnych grup elegancko połączone w jeden data frame
#' dla każdego pogrupowania, z którym związane są efekty losowe.
#' @param x model klasy \code{lmerMod}
#' @return lista data frame'ów
#' @export
#' @import lme4
ranef_ddf = function(x) {
  temp = ranef(x, condVar = TRUE)
  temp = mapply(
    function(z, nazwa) {
      postSD = matrix(NA, nrow = nrow(z), ncol = ncol(z),
                      dimnames = list(rownames(z), paste0("csd_", colnames(z))))
      for (i in 1:ncol(z)) postSD[, i] = attributes(z)$postVar[i, i, ]^0.5
      z = data.frame(rownames(z), z, postSD, check.names = FALSE)
      names(z)[1] = nazwa
      return(z)
    }, temp, as.list(names(temp)), SIMPLIFY = FALSE
  )
  return(temp)
}
