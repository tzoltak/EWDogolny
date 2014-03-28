\name{kolorki}
\alias{kolorki}
\title{Skala kolorow dla wektora liczbowego.}
\usage{
kolorki(x, kolory, zakres = NULL)
}
\arguments{
  \item{x}{wektor liczbowy}

  \item{kolory}{wektor tekstowy z nazwami kolorów}

  \item{zakres}{opcjonalny dwualementowy wektor podający
  zakres wartości, w jakim ma różnicować się skala
  kolorów}
}
\value{
wektor z nazwami kolorów o długości równej wektorowi x
}
\description{
Funkcja tworzy mapowanie wartości wektora liczbowego na
kolory (skalę kolorów) i zwraca wynik takiego mapowania.
}
\examples{
x=runif(150)
plot(x, pch=21, bg=kolorki(x, c("green", "red")))
plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red")))
plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red"), c(0.4, 0.6)))
plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red"), c(-1, 2)))
}
\seealso{
\code{\link[grDevices]{colorRampPalette}}
}

