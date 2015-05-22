#' @title Diagnostyka regresji.
#' @description
#' Funkcja rysuje wykresy diagnostyczne dla modeli regresji i w tym sensie powiela
#' funkcjonalność \code{print(lm(...))}, jednak jest zoptymalizowana do tego, żeby dawać
#' ładne wykresy dla dużych danych (i dawać je szybko).
#' Poza tym obsługuje automatyczne zapisywanie wykresów do plików PNG.
#' @details
#' Rysowane są następujące wykresy:
#' \itemize{
#'   \item{reszty w funkcji przewidywania do diagnostyki niezależności reszt od
#'   przewidywania,}
#'   \item{wykres kwantyl-kwantyl do diagnostyki normalności rozkładu reszt,}
#'   \item{pierwiastek z modułu reszt w funkcji przewidywania do diagnostyki
#'         homoscedatyczności.}
#' }
#' Dodatkowo dla modeli mieszanych efektów, dla każdego efektu losowego:
#' \itemize{
#'   \item{wykres kwantyl-kwantyl do diagnostyki normalności rozkładu BLUPsów,}
#'   \item{BLUPsy w funkcji średniej przewidywań z części stałej do diagnostyki
#'         niezależności części stałej i losowej modelu.}
#' }
#' @param model model regresji (typowo wynik działania funkcji \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}} lub \code{\link[lme4]{lmer}})
#' @param zapiszPng tekst - ścieżka do katalogu, w którym mają być zapisane pliki PNG
#' z wykresami lub NULL, jeśli wykresy mają nie być zapisywane
#' @param tytul tekst - tytuł wykresów (od tytułu będą się też rozpoczynać nazwy
#' zapisywanych plików PNG)
#' @param smoothScatter wartość logiczna - czy do rysowania wykresów rozrzutu używać
#' funkcji \code{\link[graphics]{smoothScatter}}?
#' @param loess wartość logiczna - czy do wykresów rozrzutu dodawać linię regresji
#' nieparametrycznej (estymowanej funkcją \code{\link[stats]{loess}})?
#' @param span parametr sterujący wygładzaniem regresji nieparametrycznej
#' (p. \code{\link[stats]{loess}})
#' @param control parametry regresji nieparametrycznej - lista zwrócona przez funkcję
#' \code{\link[stats]{loess.control}}
#' @param paleta funkcja definiująca paletę kolorów do wykorzystania przy rysowaniu
#' wykresów, typowo wynik wywołania funkcji \code{\link[grDevices]{colorRampPalette}}
#' @param pch jeśli smoothScatter=FALSE, znak punktu używanego do rysowania obserwacji
#' (p. \code{\link[graphics]{points}})
#' @param cex jeśli smoothScatter=FALSE, wielkość punktów używanych do rysowania obserwacji
#' (p. \code{\link[graphics]{points}})
#' @param lwd jeśli loess=TRUE, grubość rysowanej linii  (p. \code{\link[graphics]{lines}})
#' @param alpha jeśli smoothScatter=FALSE i nie odpowiada Ci wynik działania domyślnego
#' argumentu ustalajacego, jak bardzo przeźroczyste są punkty, możesz zadać przeźroczystość sam
#' @return funkcja nic nie zwraca
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[lme4]{lmer}},
#' \code{\link[graphics]{smoothScatter}}, \code{\link[stats]{loess}}
#' @examples
#' \dontrun{
#'   x = rnorm(100000)
#'   y = x + rnorm(100000)
#'   diagnostyka(lm(y ~ x), tytul = "y~x")
#'
#'   library(lme4)
#'   fm1 = lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'   diagnostyka(fm1)
#' }
#' @export
#' @import lme4
diagnostyka = function(model, zapiszPng=NULL, tytul=NULL, smoothScatter=TRUE,
                       loess=TRUE, span=0.5,
                       control = loess.control(surface = "interpolate",
                                               statistics = "approximate",
                                               trace.hat = "approximate",
                                               cell = 1),
                       paleta = colorRampPalette(c("white", blues9)),
                       pch=16, cex=1, lwd=2, alpha=NULL) {
  stopifnot(is.logical(smoothScatter), is.logical(loess),
            is.null(zapiszPng) | is.character(zapiszPng),
            is.null(tytul) | is.character(tytul),
            is.numeric(span), is.numeric(pch), is.numeric(cex), is.numeric(lwd),
            is.null(alpha) | is.numeric(alpha),
            is.function(paleta),
            is.list(control)
  )
  stopifnot(span > 0, pch > 0, cex > 0, lwd > 0, alpha > 0 & alpha <= 1)
  reszty = try(resid(model))
  przewidywania = try(fitted(model))
  if ("try-error" %in% c(class(reszty), class(przewidywania))) {
    stop("Na argumencie 'model' musi dać się wykonać funkcje resid() i fitted().")
  }

  katalog = getwd()
  if (is.character(zapiszPng)) {
    zapiszPng = try(setwd(zapiszPng[1]), silent = TRUE)
    if ("try-error" %in% class(zapiszPng)) {
      zapiszPng = NULL;
      warning("Nie udało się zapisać wykresów - podano niepoprawną nazwę (ścieżkę) katalogu.")
    }
  }
  parametryGraficzne = par(no.readonly = TRUE)
  par(bg = "white")

  if (is.null(tytul)) tytul = ""
  podtytul = paste0(paste0(deparse(formula(model), 70), collapse = "\n"))


  if (any(class(model) %in% c("lmerMod", "glmerMod"))) {
    stdReszty = resid(model, scaled = TRUE)
    blup = ranef(model)
    zagrPrzewSt = apply(t(model@pp$X) * fixef(model), 2, sum)  # chwilowo wbrew nazwie wcale nie zagregowane, ale to się zaraz zmieni
    zagrPrzewSt = mapply(  # na potrzeby diagnostyki blupsów
      function(x, y, przew) {
        return(setNames(
          aggregate(przew, list(x), mean),
          c(y, "sr_przew")
        ))
      },
      as.list(data.frame(model@frame[, names(blup)])),
      as.list(names(blup)),
      MoreArgs = list(przew = zagrPrzewSt),
      SIMPLIFY = FALSE
    )
    names(zagrPrzewSt) = names(blup)
    efLos = VarCorr(model)  # przydzadzą nam się też odchylenia standardowe efektów losowych
    efLos = lapply(efLos, function(x) return(attributes(x)$stddev))
  } else {
    stdReszty = rstandard(model)
  }
  pierwBezwzglStdReszty = sqrt(abs(stdReszty))

  kolorKropek = rgb2hsv(col2rgb(paleta(4)[3]))
  if (is.null(alpha)) alpha = ifelse(
    length(reszty) > 100,
    1/log(length(reszty)/10,10),
    1
  )
  kolorKropek = hsv(kolorKropek[1], kolorKropek[2], kolorKropek[3], alpha)

  limPrzewidywania         = range(przewidywania)
  limReszty                = max(abs(reszty)) * c(-1,1)
  limStdReszty             = max(abs(stdReszty)) * c(-1,1)
  limPierwBezwzglStdReszty = c(0, max(pierwBezwzglStdReszty))
  # reszty ~ przewidywania
  plot(NA, NA, xlim = limPrzewidywania, ylim = limReszty, main = "",
       xlab = "wartości przewidywane", ylab = "reszty")
  title(main = tytul, line = 2.2)
  title(main = list(podtytul, cex = 0.8, font = 1), line = 0.5)
  if (smoothScatter) {
    smoothScatter(przewidywania, reszty, add = TRUE, nbin = 256,
                  bandwidth = c((limPrzewidywania[2] - limPrzewidywania[1]) / 256,
                                (limReszty[2] - limReszty[1]) / 256),
                  colramp = paleta)
  } else {
    points(przewidywania, reszty, pch = pch, cex = cex, col = kolorKropek)
  }
  grid(col = grey(0.5))
  abline(h = 0)
  if (loess) {
    lines(przew_npar(reszty ~ przewidywania, data.frame(reszty, przewidywania)),
          col = 2, lwd = lwd, lty = 1)
  }
  if (is.character(zapiszPng)) {
    dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_",
                          "reszty_w_funkcji_przewidywania", ".png"),
              height = 1004, width = 1004, pointsize = 12, res = 150)
  }
  # normalność reszt
  temp = qqnorm(stdReszty, plot.it = FALSE)
  plot(NA, NA, xlim = limStdReszty, ylim = limStdReszty, main = "",
       xlab = "", ylab = "rozkład empiryczny reszt standaryzowanych")
  title(main = tytul, line = 2.2)
  title(main = list(podtytul, cex = 0.8, font = 1), line = 0.5)
  title(xlab = "rozkład teoretyczny: N(0;1)", line = 2)
  title(sub = list("poniżej przekątnej - częstość większa niż teoretyczna,\npowyżej przekątnej - częstość mniejsza niż teoretyczna",
                   cex = 0.8, font = 1), line = 3.8)
  if (smoothScatter) {
    smoothScatter(temp, add = TRUE, transformation = function(x) x^(1/8),
                  nbin = 256, bandwidth = c((limStdReszty[2] - limStdReszty[1]) / 512,
                                            (limStdReszty[2] - limStdReszty[1]) / 512),
                  colramp = paleta)
  } else {
    points(temp, pch = pch, cex = cex, col = kolorKropek)
  }
  grid(col = grey(0.5))
  abline(0, 1, col = 2, lwd = 2, lty = 1)
  if (is.character(zapiszPng)) {
    dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_",
                          "normalnosc_reszt", ".png"),
              height = 1004, width = 1004, pointsize = 12, res = 150)
  }
  # homoscedatyczność
  plot(NA, NA, xlim = limPrzewidywania, ylim = limPierwBezwzglStdReszty,
       main = "", xlab = "wartości przewidywane", ylab = "")
  title(main = tytul, line = 2.2)
  title(main = list(podtytul, cex = 0.8, font = 1), line = 0.5)
  title(ylab = list(as.expression(substitute(sqrt(abs(R)),
                                             list(R = as.name("reszty standaryzowane"))))),
        line = 2.5)
  if (smoothScatter) {
    smoothScatter(przewidywania, pierwBezwzglStdReszty, add = TRUE, nbin = 256,
                  bandwidth = c((limPrzewidywania[2] - limPrzewidywania[1]) / 256,
                                (limPierwBezwzglStdReszty[2] - limPierwBezwzglStdReszty[1]) / 256),
                  colramp = paleta)
  } else {
    points(przewidywania, pierwBezwzglStdReszty, pch = pch, cex = cex, col = kolorKropek)
  }
  grid(col = grey(0.5))
  if (loess) {
    lines(przew_npar(pierwBezwzglStdReszty ~ przewidywania,
                     data.frame(pierwBezwzglStdReszty, przewidywania)),
          col = 2, lwd = lwd, lty = 1)
  }
  if (is.character(zapiszPng)) {
    dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_",
                          "homoscedatycznosc", ".png"),
              height = 1004, width = 1004, pointsize = 12, res = 150)
  }
  # diagnostyka efektów losowych (z poziomów innych niż indywidualny)
  if (any(class(model) %in% c("lmerMod", "glmerMod"))) {
    for (i in 1:length(blup)) {
      for (j in 1:ncol(blup[[i]])) {
        blupTemp = blup[[i]][[j]]
        zagrPrzewStTemp = zagrPrzewSt[names(blup)[[i]]][[1]]$sr_przew
        # analiza normalności rozkładów BLUPsów;
        skala = unlist(lapply(efLos[grepl(paste0("^", names(blup)[i],
                                                 "(|[.][[:digit:]]+)$"),
                                          names(efLos))],
                              function(x, z) {return(x[names(x) == z])},
                              z = names(blup[[i]])[j]))
        temp = qqnorm(blupTemp / skala, plot.it = FALSE)
        limBlup = max(abs(unlist(temp))) * c(-1,1)
        plot(NA, NA, xlim = limBlup, ylim = limBlup, main = "", xlab = "",
             ylab = paste0("rozkł. emp. efektu losowego: ",
                           names(blup[[i]])[j], " | ", names(blup)[i]))
        title(main = tytul, line = 2.2)
        title(main = list(podtytul, cex = 0.8, font = 1), line = 0.5)
        title(xlab = "rozkład teoretyczny: N(0;1)", line = 2)
        title(sub = list("poniżej przekątnej - częstość większa niż teoretyczna,\npowyżej przekątnej - częstość mniejsza niż teoretyczna",
                         cex = 0.8, font = 1), line = 3.8)
        if (smoothScatter) {
          smoothScatter(temp, add = TRUE, transformation = function(x) x^(1/8),
                        nbin = 256, bandwidth = c((limBlup[2] - limBlup[1]) / 512,
                                                (limBlup[2] - limBlup[1]) / 512),
                        colramp = paleta)
        } else {
          points(temp, pch = pch, cex = cex, col = kolorKropek)
        }
        grid(col = grey(0.5))
        abline(0, 1, col = 2, lwd = 2, lty = 1)
        if (is.character(zapiszPng)) {
          dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)),
                                "_", "normalnosc_ef_los_", i, "_", j, ".png"),
                    height = 1004, width = 1004, pointsize = 12, res = 150)
        }
        # analiza związków między średnimi przewidywanami wynikającymi z części stałej modelu a BLUPsami;
        limBlup = max(abs(unlist(blupTemp))) * c(-1,1)
        limPrzewidywaniaZagr = range(zagrPrzewStTemp)
        if (limPrzewidywaniaZagr[2] > limPrzewidywaniaZagr[1]) {  # dla wyników badań eksperymentalnych ze ślicznie zrównoważonymi schematami może być tak, że nie ma tu czego sprawdzać
          plot(NA, NA, xlim = limPrzewidywaniaZagr, ylim = limBlup, main = "",
               xlab = paste0("śr. wartości przewidywane w ramach grup\nwyróżnionych ze względu na zmienną: ",
                             names(blup)[i]),
               ylab = paste0("przew. dla efektu losowego: ",
                             names(blup[[i]])[j], " | ", names(blup)[i]))
          title(main = tytul, line = 2.2)
          title(main = list(podtytul, cex = 0.8, font = 1), line = 0.5)
          if (smoothScatter) {
            smoothScatter(zagrPrzewStTemp, blupTemp, add = TRUE, nbin = 256,
                          bandwidth = c((limPrzewidywaniaZagr[2] - limPrzewidywaniaZagr[1]) / 256,
                                        (limBlup[2] - limBlup[1]) / 256),
                          colramp = paleta)
          } else {
            points(zagrPrzewStTemp, blupTemp, pch = pch, cex = cex, col = kolorKropek)
          }
          grid(col = grey(0.5))
          abline(h = 0)
          if (loess) {
            lines(przew_npar(blupTemp ~ zagrPrzewStTemp,
                             data.frame(blupTemp, zagrPrzewStTemp)),
                  col = 2, lwd = lwd, lty = 1)
          }
          if (is.character(zapiszPng)) {
            dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)),
                                  "_", "przew_ef_los_w_funkcji_przewidywania",
                                  i, "_", j, ".png"),
                      height = 1004, width = 1004, pointsize = 12, res = 150)
          }
        } else {
          message(paste0("Dla pogrupowania ze względu na zmienną '",
                         names(blup)[i], "' nie wygenerowano wykresu\n   przewidywania ef. los. | przewidywanie z ef. stałych\nze względu na dosknale zrównoważony schemat analizy (i w efekcie brak zróżnicowania przewidywania z ef. stałych w ramach grup)."))
        }
      }
    }
  }
  # kończenie
  par(parametryGraficzne)
  setwd(katalog)
  invisible(NULL)
}
