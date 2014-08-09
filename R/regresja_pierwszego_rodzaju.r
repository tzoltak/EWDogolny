#' @title Regresja pierwszego rodzaju.
#' @description
#' Funkcja wylicza regresję pierwszego rodzaju (inaczej mówiąc, parametry warunkowe), domyślnie średnich.
#' @param zmZal ciąg znaków - nazwa zmiennej zależnej
#' @param zmNZal wektor tekstowy z nazwami zmiennych niezależnych
#' @param dane ramka danych zawierająca zmienne występujące w formule
#' @param FUN funkcja służąca do wyliczeniu parametrów warunkowych
#' @details
#' Funkcja \code{FUN} powinna zwracać pojedynczą wartość liczbową i przyjmować argument \code{na.rm}.
#' @return data frame
#' @import plyr
#' @export
regr_pierw_rodz = function(zmZal, zmNZal, dane, FUN=mean) {
	stopifnot(is.data.frame(dane), is.function(FUN),
						is.character(zmZal ), length(zmZal ) == 1,
						is.character(zmNZal), length(zmNZal) >  0
	)
	stopifnot(zmZal %in% names(dane), all(zmNZal %in% names(dane)))
	temp = ddply(dane, zmNZal,
							 function(x, y, FUN) {
							 	return(do.call(FUN, list(x=x[, y], na.rm=TRUE)))
							 },
							 y=zmZal, FUN=FUN
	)
	names(temp)[ncol(temp)] = zmZal
	return(temp)
}
