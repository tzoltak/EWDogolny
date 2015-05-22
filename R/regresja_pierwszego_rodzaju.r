#' @title Regresja pierwszego rodzaju.
#' @description
#' Funkcja wylicza regresję pierwszego rodzaju (inaczej mówiąc, parametry
#' warunkowe), domyślnie średnich.
#' @param zmZal ciąg znaków - nazwa zmiennej zależnej
#' @param zmNZal wektor tekstowy z nazwami zmiennych niezależnych
#' @param dane ramka danych zawierająca zmienne występujące w formule
#' @param FUN funkcja służąca do wyliczeniu parametrów warunkowych
#' @details
#' Funkcja \code{FUN} powinna zwracać pojedynczą wartość liczbową i przyjmować
#' argument \code{na.rm}.
#' @return data frame
#' @import plyr
#' @export
regr_pierw_rodz = function(zmZal, zmNZal, dane, FUN = mean) {
	stopifnot(is.data.frame(dane), is.function(FUN),
						is.character(zmZal ), length(zmZal ) == 1,
						is.character(zmNZal), length(zmNZal) >  0
	)
	stopifnot(zmZal %in% names(dane), all(zmNZal %in% names(dane)))
	temp = ddply(dane, zmNZal,
							 function(x, y, FUN) {
							 	return(do.call(FUN, list(x = x[, y], na.rm = TRUE)))
							 },
							 y = zmZal, FUN = FUN
	)
	names(temp)[ncol(temp)] = zmZal
	return(temp)
}
#' @title Regresja nieparametryczna.
#' @description
#' Funkcja przy pomocy funkcji \code{\link[stats]{loess}} wylicza dopasowanie
#' nieparametryczne i przetwarza wyniki w formę zdatnę do łatwego rysowania.
#' @details
#' Funkcja wywołuje \code{\link[stats]{loess}} z parametrami ustawionymi tak,
#' aby nie liczyło się to za długo. Przy tym troszeczkę dopasowuje dokładność
#' (konkretnie parametr \code{cell} do liczby obserwacji, żeby nie przesadzić
#' z upraszczaniem.)
#' @param formula formuła z modelem zależności
#' @param data ramka danych zawierająca zmienne występujące w formule
#' @param degree parametr sterujący wygładzaniem (p. \code{\link[stats]{loess}})
#' @param span parametr sterujący wygładzaniem (p. \code{\link[stats]{loess}})
#' @param control lista zwrócona przez funkcję \code{\link[stats]{loess.control}}
#' @return ramka danych zawierająca zmienne 'x' (wartości zm. niezależnej) i 'y'
#' (przewidywanie)
#' @seealso \code{\link[stats]{loess}}, \code{\link[stats]{loess.control}}
#' @examples
#' x=runif(100000) * 4 - 2
#' y=x+rnorm(100000)
#' z=x+rnorm(100000)
#' plot(z, y, pch='.', cex=0.5)
#' grid(col=grey(0.5))
#' lines(przew_npar(y~z, data.frame(z, y)), col=2, lwd=2)
#' @export
przew_npar = function(formula, data, degree = 2, span = 0.5, control = NULL) {
  stopifnot(
    "formula" %in% class(formula), is.data.frame(data) | is.matrix(data),
    degree %in% c(0,1,2), is.numeric(span), is.null(control) | is.list(control)
  )
  stopifnot(span > 0, all( all.vars(formula) %in% colnames(data) ))

  # w zależności od liczby obserwacji dopasowujemy prametr 'cell'
  # aby dla dużych nie liczyło się za długo, a dla krótkich nie liczyło się zbyt zgrubnie
  if (is.null(control)) {
    if (nrow(data) > 100000) {
      cell = 1
    } else if (nrow(data) > 10000) {
      cell = 0.5
    } else {
      cell = 0.2
    }
    control = loess.control(surface = "interpolate", statistics = "approximate",
                            trace.hat = "approximate", cell = cell)
  }
  # sama estymacja i obróbka wyniku
  przewNpar = loess(formula, data, degree = degree, span = span, control = control)
  przewNpar = unique(setNames(
    as.data.frame(przewNpar[c("x", "fitted")]),
    c("x", "y")
  ))
  return(przewNpar[order(przewNpar$x), ])
}
#' @title Regresja nieparametryczna.
#' @description
#' Funkcja powiela działanie funkcji \code{\link{przew_npar}} i jest uruchamiana
#' z parametrami funkcji \code{\link{regr_pierw_rodz}}.
#' @param zmZal ciąg znaków - nazwa zmiennej zależnej
#' @param zmNZal wektor tekstowy z nazwami zmiennych niezależnych
#' @param dane ramka danych zawierająca zmienne występujące w formule
#' @param przew_nparPar lista z dodatkowymi parametrami dla funkcji przew_npar.
#' @return ramka danych zawierająca zmienne 'x' (wartości zm. niezależnej) i 'y'
#' (przewidywanie)
#' @export
przew_npar_rpr = function(zmZal, zmNZal, dane, przew_nparPar = NULL){
  if (is.null(przew_nparPar)){
    przew_nparPar = list()
  }
  przew_nparPar$formula = as.formula(paste0(zmZal, "~", zmNZal))
  przew_nparPar$data = dane
  przewNPar = do.call(przew_npar, przew_nparPar)
  return(przewNPar)
}
