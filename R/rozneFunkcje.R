#' @title Skala kolorow dla wektora liczbowego.
#' @description
#' Funkcja tworzy mapowanie wartości wektora liczbowego na kolory (skalę kolorów)
#' i zwraca wynik takiego mapowania.
#' @param x wektor liczbowy
#' @param kolory wektor tekstowy z nazwami kolorów
#' @param zakres opcjonalny dwualementowy wektor podający zakres wartości,
#' w jakim ma różnicować się skala kolorów
#' @return wektor z nazwami kolorów o długości równej wektorowi x
#' @seealso \code{\link[grDevices]{colorRampPalette}}
#' @examples
#' x=runif(150)
#' plot(x, pch=21, bg=kolorki(x, c("green", "red")))
#' plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red")))
#' plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red"), c(0.4, 0.6)))
#' plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red"), c(-1, 2)))
#' @export
kolorki = function(x, kolory, zakres=NULL) {
  stopifnot(is.numeric(x), is.character(kolory), is.null(zakres) | is.numeric(zakres))
  stopifnot(length(kolory) > 1, is.null(zakres) | length(zakres) == 2)

  if (is.null(zakres)) {  # potnij x na 100 przedziałów o równej szerokości
    x = cut(x, c(-Inf,
                 seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
                     length.out = 101)[-c(1,101)],
                 Inf))
  } else {
    # potnij x wg przedziałów powstałych przez podział podanego zakresu
    # na 100 kawałków o równej szerokości
    x = cut(x, c(-Inf, seq(zakres[1], zakres[2], length.out = 101)[-c(1,101)], Inf))
  }
  paleta = colorRampPalette(kolory)
  lKolorow = length(unique(x))
  return(paleta(nlevels(x))[x])
}
#' @title Konwersja nazw OKE
#' @description
#' Funkcja konwertuje nazwy OKE pomiędzy dwoma standardami stosowanymi w zbiorach EWD:
#' 1) zapisem z użyciem polskich znaków i wielkich liter,
#' 2) zapisem tylko z wykorzystaniem małych liter, bez polskich znaków.
#' @param x czynnik opisujący przydział do OKE
#' @param naPolskie wartość logiczna - czy konwertować do standardu zapisu
#' z polskimi znakami? (jeśli FALSE - konwersja do standardu zapisu bez polskich znaków)
#' @return czynnik
#' @examples
#' x=c("gdansk", "lodz")
#' y=c("Warszawa", "Kraków")
#' nazwy_oke(factor(x))
#' nazwy_oke(factor(y), FALSE)
#' nazwy_oke(factor(c(x, y)))
#' @export
nazwy_oke = function(x, naPolskie=TRUE) {
  stopifnot(is.factor(x), is.logical(naPolskie))
  stopifnot(length(naPolskie) == 1, naPolskie %in% c(TRUE, FALSE))
  stopifnot(all(levels(x) %in% c("gdansk", "jaworzno", "krakow", "lomza", "lodz",
                                 "poznan", "warszawa", "wroclaw", "Gdańsk",
                                 "Jaworzno", "Kraków", "Łomża", "Łódź", "Poznań",
                                 "Warszawa", "Wrocław")))
  if (!(all(levels(x) %in% c("gdansk", "jaworzno", "krakow", "lomza", "lodz",
                             "poznan", "warszawa", "wroclaw")) |
        all(levels(x) %in% c("Gdańsk", "Jaworzno", "Kraków", "Łomża", "Łódź",
                             "Poznań", "Warszawa", "Wrocław")))) {
    warning("Wykryto poziomy zarówno w formacie z polskimi znakami jak i w formacie bez polskich znaków.")
  }
  if (naPolskie) {
    levels(x) = sub("gdansk",   "Gdańsk",   levels(x))
    levels(x) = sub("jaworzno", "Jaworzno", levels(x))
    levels(x) = sub("krakow",   "Kraków",   levels(x))
    levels(x) = sub("lomza",    "Łomża",    levels(x))
    levels(x) = sub("lodz",     "Łódź",     levels(x))
    levels(x) = sub("poznan",   "Poznań",   levels(x))
    levels(x) = sub("warszawa", "Warszawa", levels(x))
    levels(x) = sub("wroclaw",  "Wrocław",  levels(x))
    x = factor(as.character(levels(x))[x],
               levels = c("Gdańsk", "Jaworzno", "Kraków", "Łomża", "Łódź",
                          "Poznań", "Warszawa", "Wrocław"))
  } else {
    levels(x) = sub("Gdańsk",   "gdansk",   levels(x))
    levels(x) = sub("Jaworzno", "jaworzno", levels(x))
    levels(x) = sub("Kraków",   "krakow",   levels(x))
    levels(x) = sub("Łomża",    "lomza",    levels(x))
    levels(x) = sub("Łódź",     "lodz",     levels(x))
    levels(x) = sub("Poznań",   "poznan",   levels(x))
    levels(x) = sub("Warszawa", "warszawa", levels(x))
    levels(x) = sub("Wrocław",  "wroclaw",  levels(x))
    x = factor(as.character(levels(x))[x],
               levels = c("gdansk", "jaworzno", "krakow", "lomza", "lodz",
                          "poznan", "warszawa", "wroclaw"))
  }
  return(x)
}
#' @title Wazone kwantyle
#' @description
#' Funkcja wylicza kwantyle uwzględniając nierówne wagi jednostek obserwacji.
#' Sposób wyliczania kwanrtyli odpowiada podejściu \code{type=2} z funkcji
#' \code{\link[stats]{quantile}}.
#' @param x wektor wartości zmiennej
#' @param wagi wektor wag
#' @param probs wektor prawdopodobieństw definiujących kwantyle, które mają
#' zostać wyliczone
#' @param na.rm wartość logiczna - czy pominąć braki danych?
#' @return wektor z wyliczonymi kwantylami
#' @examples
#' x=1:10
#' w=rep(1:2, each=5)
#' kwantyl_wazony(x, w)
#' quantile(c(1:5, rep(6:10, each=2)), type=2)
#' @export
kwantyl_wazony = function(x, wagi, probs=seq(0, 1, 0.25), na.rm=FALSE) {
  stopifnot(is.numeric(x), is.numeric(wagi), is.numeric(probs), is.logical(na.rm))
  stopifnot(all(wagi >= 0), length(wagi) == length(x), all(probs >= 0),
            all(probs <= 1), length(na.rm) == 1)

  temp = data.frame(x = x, wagi = wagi, wagiOdwr = wagi)
  # obsługa ew. braków danych w x
  if (na.rm) {  # gdzie braki, to zignoruj i po prostu policz na pozostałych
    temp = na.omit(temp)
  } else {
    if (with(temp, any(is.na(x[wagi > 0])))) {  # jeśli jakieś obserwacje z niezerowymi wagami mają braki, zwróć brak danych
      return(setNames(rep(NA, length(probs)), probs))
    } else {  # jeśli braki dotyczą tylko obserwacji z wagami równymi zero, można je zignorować
      temp = na.omit(temp)
    }
  }
  # samo wyliczanie kwantyli
  temp = temp[order(temp$x), ]
  temp$wagi[!is.na(temp$wagi)] =
    cumsum(temp$wagi[!is.na(temp$wagi)] ) / sum(temp$wagi, na.rm = TRUE)
  temp$wagiOdwr[!is.na(temp$wagiOdwr)] =
    rev(cumsum(rev(temp$wagiOdwr[!is.na(temp$wagiOdwr)]))) / sum(temp$wagiOdwr, na.rm = TRUE)
  wynik = vector()
  for (i in probs) {
    wynik = c(wynik, mean(temp$x[temp$wagi >= i & temp$wagiOdwr >= (1 - i)]))
  }
  return(setNames(wynik, probs))
}
