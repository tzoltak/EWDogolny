% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regresja_pierwszego_rodzaju.R
\name{przew_npar}
\alias{przew_npar}
\title{Regresja nieparametryczna.}
\usage{
przew_npar(formula, data, degree = 2, span = 0.5, control = NULL)
}
\arguments{
\item{formula}{formuła z modelem zależności}

\item{data}{ramka danych zawierająca zmienne występujące w formule}

\item{degree}{parametr sterujący wygładzaniem (p. \code{\link[stats]{loess}})}

\item{span}{parametr sterujący wygładzaniem (p. \code{\link[stats]{loess}})}

\item{control}{lista zwrócona przez funkcję \code{\link[stats]{loess.control}}}
}
\value{
ramka danych zawierająca zmienne 'x' (wartości zm. niezależnej) i 'y'
(przewidywanie)
}
\description{
Funkcja przy pomocy funkcji \code{\link[stats]{loess}} wylicza dopasowanie
nieparametryczne i przetwarza wyniki w formę zdatnę do łatwego rysowania.
}
\details{
Funkcja wywołuje \code{\link[stats]{loess}} z parametrami ustawionymi tak,
aby nie liczyło się to za długo. Przy tym troszeczkę dopasowuje dokładność
(konkretnie parametr \code{cell} do liczby obserwacji, żeby nie przesadzić
z upraszczaniem.)
}
\examples{
x=runif(100000) * 4 - 2
y=x+rnorm(100000)
z=x+rnorm(100000)
plot(z, y, pch='.', cex=0.5)
grid(col=grey(0.5))
lines(przew_npar(y~z, data.frame(z, y)), col=2, lwd=2)
}
\seealso{
\code{\link[stats]{loess}}, \code{\link[stats]{loess.control}}
}

