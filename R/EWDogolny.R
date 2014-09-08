#####################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!!!! Jeśli coś tu zmieniasz, nie zapomnij sprawdzić,       !!!!! #
# !!!!! czy po modyfikacji wszystko działa.                   !!!!! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
#####################################################################

#' @title Eksport danych do schowka systemowego.
#' @description
#' Funkcja przetwarza listę składającą się z ramek danych, macierzy,i/lub wektorów na jedną wielką macierz i kopiuje ją do schowka systemowego, skąd można wkleić ją w inne miejsca.
#' Stara się przy tym dosyć inteligentnie obchodzić z nazwami (elemetntów listy, wierszy, kolumn, elementów wektorów), chociaż pewnie nie zawsze jej się udaje.
#' @param x obiekt, który ma być wyeksportowany
#' @param nw czy eksportować nazwy wierszy macierzy?
#' @param oel wartość logiczna - czy na końcu całości dodać przeniesienie do nowej linii? (nie ruszać bez potrzeby)
#' @param sep separator pola przy zapisie (jak w \code{\link[utils]{write.table}}; nie ruszać bez potrzeby)
#' @param dec separator dziesiętny przy zapisise (jak w \code{\link[utils]{write.table}}; nie ruszać bez potrzeby)
#' @param na ciąg znaków na oznaczenie braków danych przy zapisie (jak w \code{\link[utils]{write.table}}; nie ruszać bez potrzeby)
#' @param plik nazwa strumienia, do którego ma pójść zapis (nie ruszać bez potrzeby)
#' @return funkcja nic nie zwraca
#' @seealso \code{\link[base]{file}}, \code{\link[utils]{write.table}}
#' @examples
#' x = c(1:10)
#' y = 2 * x + 5
#' z = letters[11:20]
#' lista = list("ramka danych"=data.frame(x, y, z), macierz=matrix(c(x, y), ncol=2, dimnames=list(z, c("a", "b"))), wektor=setNames(x, z))
#' do_schowka(lista)
#' do_schowka(lista, nw=TRUE)
#' @export
do_schowka = function(x, nw=FALSE, oel=TRUE, sep="\t", dec=",", na="", plik="clipboard-128") {
	temp = file(plik, "w")
	options(warn=-1)
	try(cat("", file=temp))
	if (!is.list(x) | is.data.frame(x)) x = list(x)
	for (i in 1:length(x)) {
		if (!is.null(names(x))) if (names(x)[i] != "") try(write(names(x)[i], temp, append=TRUE))
		if (nw & (is.data.frame(x[[i]]) | is.matrix(x[[i]]))) {
			x[[i]] = data.frame(" "=rownames(x[[i]]), x[[i]], check.names=FALSE)
		} else if (is.matrix(x[[i]])) {
			x[[i]] = as.data.frame(x[[i]])
		}
		if (!is.null(names(x[[i]]))) try(write(names(x[[i]]), temp, append=TRUE, ncolumns=length(names(x[[i]])), sep=sep))
		if (is.vector(x[[i]])) x[[i]] = as.data.frame(as.list(x[[i]]))
		try(write.table(x[[i]], temp, sep=sep, dec=dec, na=na, row.names=FALSE, col.names=FALSE, append=TRUE))
		if (oel & i < length(x)) try(cat("\n", file=temp, append=TRUE))
	}
	close(temp)
	options(warn=0)
	invisible(NULL)
}
#' @title Podsumowanie informacji o modelach regresjii liniowej.
#' @description
#' Funkcja przygotowuję listę z parametrami modelu (modeli) regresji MNK wyestymowanych funkcją \code{\link[stats]{lm}} (na innych nie gwarantuje się poprawności działania funkcji).
#' @details
#' Dla każdego modelu zwracane są trzy tabele:
#' \itemize{
#' \item zawierająca wskaźniki dopasowania (R2, skoryg. R2, parametry testu F, AIC, BIC),
#' \item zawierająca parametry modelu, w tym standaryzowane ("xy" - normalna standaryzacja i "x" - coś na podobieństwo d Cohena),
#' \item zawierająca macierz korelacji efektów.
#' }
#' @param modele model lub lista modeli
#' @param pokaz wartość logiczna - czy wyświetlić informacje na konsoli?
#' @return lista
#' @seealso \code{\link[lme4]{lmer}}, \code{\link{do_schowka}}
#' @examples
#' x=rnorm(100)
#' y=x+rnorm(100)
#' z=x+rnorm(100,0,2)
#' modele=list("model y"=lm(y~x), "model z"=lm(z~x))
#' do_schowka(przygotuj_do_schowka_es(modele))
#' @export
przygotuj_do_schowka_es = function(modele, pokaz=TRUE) {
	stopifnot(is.logical(pokaz))
	if (!is.list(modele) | ("lm" %in% class(modele))) modele = list(modele)
	return(unlist(lapply(modele, 
		function(x) {
			y = try(summary(x, correlation=TRUE))
			if (pokaz) {cat("\n#############################################"); print(y)}
			if ("summary.lm" %in% class(y)) {
				# macierz danych można wyciągnąć odwracając dekompozycję QR
				daneM = data.frame(x$model[,1], qr.X(x$qr)[, -1], check.names=FALSE)
				odchStd = unlist(lapply(daneM, sd))
				return(list(
					dopasowanie = c(
						R2           = y$r.squared,
						"skoryg. R2" = y$adj.r.squared,
						F            = setNames(y$fstatistic[1], ""),
						"F df1"      = setNames(y$fstatistic[2], ""),
						"F df2"      = setNames(y$fstatistic[3], ""),
						AIC          = AIC(x),
						BIC          = BIC(x)
					),
					efSt = data.frame(
						param             = rownames(y$coefficients), y$coefficients,
						"std.xy Estimate" = y$coefficients[, 1] * c(NA, odchStd[2:nrow(y$coefficients)]) / odchStd[1],
						"std.x Estimate"  = y$coefficients[, 1] / odchStd[1],
						check.names=FALSE
					),
					# wstawka z do.call w celu wyczyszczenia przekątnej i górnego trójkąta macierzy korelacji
					corEfSt = do.call(
						function(x) {
							x[!lower.tri(x)] = NA
							return(data.frame(" "=rownames(x), x, check.names=FALSE))
						},
						list(y$correlation)
					)
				))
			} else {
				warning("Obiekt nie jest klasy 'lm'.")
				return(NULL)
			}
		}
	), recursive=FALSE))
}
#' @title Podsumowanie informacji o modelach regresjii mieszanych efektow.
#' @description
#' Funkcja przygotowuję listę z parametrami modelu (modeli) regresji mieszanych efektów wyestymowanych funkcją \code{\link[lme4]{lmer}}.
#' @details
#' Dla każdego modelu zwracane są cztery tabele:
#' \itemize{
#' \item zawierająca wskaźniki dopasowania,
#' \item zawierającą macierz wariancji/kowariancji efektów losowych,
#' \item zawierająca parametry modelu, w tym standaryzowane ("xy" - normalna standaryzacja i "x" - coś na podobieństwo d Cohena),
#' \item zawierająca macierz korelacji efektów stałych.
#' }
#' Standaryzacja dokonywana jest względem wariancji związanej z częścią stałą modelu i błędami indywidualnymi.
#' @param modele model lub lista modeli
#' @param pokaz wartość logiczna - czy wyświetlić informacje na konsoli?
#' @return lista
#' @seealso \code{\link[stats]{lm}}, \code{\link{do_schowka}}
#' @examples
#' library(lme4)
#' fm1 = lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' fm2 = lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
#' do_schowka(przygotuj_do_schowka_me(list(fm1, fm2)))
#' @export
przygotuj_do_schowka_me = function(modele, pokaz=TRUE) {
	stopifnot(is.logical(pokaz))
	if (!is.list(modele)) modele = list(modele)
	return(unlist(lapply(modele, 
		function(x) {
			y = summary(x)
			if (pokaz) {
				cat("\n#############################################\n")
				print(y, correlation=TRUE)
			}
			if ("lmerMod" %in% class(x)) {
				# w obiektami klasy mer lub summary.mer jest prościej, niż z lm
				odchStd = unlist(lapply(as.data.frame(x@pp$X[, -1]), sd))
				# standaryzacja idzie względem wariancji zm. zależnej na poziomie indywidualnym (wariancja przewidywania + wariancja błędów indywidualnych)
				warPrzew = var(apply(t(x@pp$X[, -1]) * fixef(x)[-1], 2, sum), na.rm=TRUE)
				odchStd = c(zal=sqrt(warPrzew + attributes(VarCorr(x))$sc^2), odchStd)
				# za to z efektami losowymi trzeba się napocić
				efLos = mapply(
					function(z, nazwa) {
						return(data.frame(
							"zm. grupująca" = c(nazwa, rep("", nrow(z) - 1)),
							"efekt dla"     = colnames(z),
							"wariancja"     = attributes(z)$stddev^2,
							"odch. stand."  = attributes(z)$stddev,
							do.call(
								function(x) {
									x[!lower.tri(x)] = NA
									colnames(x)=c("cor", rep("", ncol(x)-1))
									return(x)
								},
								list(x=attributes(z)$correlation)
							),
							check.names=FALSE
						)[, - (4 + nrow(z))])
					},
					VarCorr(x),
					as.list(names(VarCorr(x))),
					SIMPLIFY=FALSE
				)
				maxLK = max(unlist(lapply(efLos, ncol)))
				efLos = lapply(efLos,
					function(z, maxLK) {
						if (ncol(z) < maxLK) z = cbind(z, matrix(NA, nrow=nrow(z), ncol=maxLK - ncol(z)))
						if (ncol(z) > 4) names(z)[5:ncol(z)] = c("cor", rep("", ncol(z)-5))
						return(z)
					}, maxLK=maxLK
				)
				if (length(efLos) > 1) {
					temp=rbind(efLos[[1]], efLos[[2]])
				} else if (length(efLos) == 1) {
					temp = efLos[[1]]
				}
				if (length(efLos) >  2) {
					for (i in 3:length(efLos)) temp = rbind(temp, efLos[[i]])
				}
				if (length(efLos) == 0) {
					efLos = NULL
				} else {
					efLos=rbind(
						temp,
						setNames(
							data.frame(
								c("reszty ind.", "ef. stałe"),
								"",
								c(attributes(VarCorr(x))$sc^2, warPrzew),
								c(attributes(VarCorr(x))$sc  , warPrzew^0.5),
								matrix(NA, nrow=2, ncol=ncol(temp)-4)
							),
							names(temp))
					)
					rownames(efLos) = NULL
				}
				# końcówka				
				return(list(
					AICtab = data.frame(
						kryterium = attributes(y$AICtab)$names,
						"l. obs." = attributes(logLik(x))$nobs,
						LL        = setNames(logLik(x), ""),
						df        = attributes(logLik(x))$df,
						deviance  = setNames(deviance(x),""),
						check.names=FALSE
					),
					efLos  = efLos,
					efSt   = data.frame(
						param             = rownames(coef(y)), coef(y),
						"std.xy Estimate" = coef(y)[, 1] * c(NA, odchStd[2:nrow(coef(y))]) / odchStd[1],
						"std.x Estimate"  = coef(y)[, 1] / odchStd[1],
						check.names=FALSE),
					# wstawka z do.call w celu wyczyszczenia przekątnej i górnego trójkąta macierzy korelacji
					corEfSt=do.call(
						function(x) {
							x[!lower.tri(x)]=NA
							return(data.frame(
								" " = rownames(x),
								x,
								check.names=FALSE
							))
						},
						list(matrix(
							vcov(x)@factors$correlation@x,
							nrow=vcov(x)@factors$correlation@Dim[1],
							dimnames=list(
								rownames(coef(y)),
								rownames(coef(y))
							)
						))
					)
				))
			} else {
				warning("Obiekt nie jest klasy 'lmerMod'.")
				return(NULL)
			}
		}
	), recursive=FALSE))
}
#' @title Regresja nieparametryczna.
#' @description
#' Funkcja powiela działanie funkcji \code{\link{przew_npar}} i jest uruchamiana z parametrami funkcji \code{\link{regr_pierw_rodz}}.
#' @param zmZal ciąg znaków - nazwa zmiennej zależnej
#' @param zmNZal wektor tekstowy z nazwami zmiennych niezależnych
#' @param dane ramka danych zawierająca zmienne występujące w formule
#' @param przew_nparPar lista z dodatkowymi parametrami dla funkcji przew_npar.
#' @return ramka danych zawierająca zmienne 'x' (wartości zm. niezależnej) i 'y' (przewidywanie)
#' @export
przew_npar_rpr = function(zmZal, zmNZal, dane, przew_nparPar = NULL){
  if(is.null(przew_nparPar)){
    przew_nparPar = list()
  }
  przew_nparPar$formula = as.formula(paste0(zmZal, "~", zmNZal))
  przew_nparPar$data = dane
  przewNPar = do.call(przew_npar, przew_nparPar)
}
#' @title Regresja nieparametryczna.
#' @description
#' Funkcja przy pomocy funkcji \code{\link[stats]{loess}} wylicza dopasowanie nieparametryczne i przetwarza wyniki w formę zdatnę do łatwego rysowania.
#' @details
#' Funkcja wywołuje \code{\link[stats]{loess}} z parametrami ustawionymi tak, aby nie liczyło się to za długo. Przy tym troszeczkę dopasowuje dokładność (konkretnie parametr \code{cell} do liczby obserwacji, żeby nie przesadzić z upraszczaniem.)
#' @param formula formuła z modelem zależności
#' @param data ramka danych zawierająca zmienne występujące w formule
#' @param degree parametr sterujący wygładzaniem (p. \code{\link[stats]{loess}})
#' @param span parametr sterujący wygładzaniem (p. \code{\link[stats]{loess}})
#' @param control lista zwrócona przez funkcję \code{\link[stats]{loess.control}}
#' @return ramka danych zawierająca zmienne 'x' (wartości zm. niezależnej) i 'y' (przewidywanie)
#' @seealso \code{\link[stats]{loess}}, \code{\link[stats]{loess.control}}
#' @examples
#' x=runif(100000) * 4 - 2
#' y=x+rnorm(100000)
#' z=x+rnorm(100000)
#' plot(z, y, pch='.', cex=0.5)
#' grid(col=grey(0.5))
#' lines(przew_npar(y~z, data.frame(z, y)), col=2, lwd=2)
#' @export
przew_npar = function(formula, data, degree=2, span=0.5, control=NULL) {
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
		control=loess.control(surface="interpolate", statistics="approximate", trace.hat="approximate", cell=cell)
	}
	# sama estymacja i obróbka wyniku
	przewNpar = loess(formula, data, degree=degree, span=span, control=control)
	przewNpar = unique(setNames(
		as.data.frame(przewNpar[c("x", "fitted")]),
		c("x", "y")
	))
	return(przewNpar[order(przewNpar$x), ])
}
#' @title Diagnostyka regresji.
#' @description
#' Funkcja rysuje wykresy diagnostyczne dla modeli regresji i w tym sensie powiela funkcjonalność \code{print(lm(...))}, jednak jest zoptymalizowana do tego, żeby dawać ładne wykresy dla dużych danych (i dawać je szybko).
#' Poza tym obsługuje automatyczne zapisywanie wykresów do plików PNG.
#' @details
#' Rysowane są następujące wykresy:
#' \itemize{
#' \item reszty w funkcji przewidywania do diagnostyki niezależności reszt od przewidywania,
#' \item wykres kwantyl-kwantyl do diagnostyki normalności rozkładu reszt,
#' \item pierwiastek z modułu reszt w funkcji przewidywania do diagnostyki homoscedatyczności.
#' }
#' Dodatkowo dla modeli mieszanych efektów, dla każdego efektu losowego:
#' \itemize{
#' \item wykres kwantyl-kwantyl do diagnostyki normalności rozkładu BLUPsów,
#' \item BLUPsy w funkcji średniej przewidywań z części stałej do diagnostyki niezależności części stałej i losowej modelu.
#' }
#' @param model model regresji (typowo wynik działania funkcji \code{\link[stats]{lm}}, \code{\link[stats]{glm}} lub \code{\link[lme4]{lmer}})
#' @param zapiszPng tekst - ścieżka do katalogu, w którym mają być zapisane pliki PNG z wykresami lub NULL, jeśli wykresy mają nie być zapisywane
#' @param tytul tekst - tytuł wykresów (od tytułu będą się też rozpoczynać nazwy zapisywanych plików PNG)
#' @param smoothScatter wartość logiczna - czy do rysowania wykresów rozrzutu używać funkcji \code{\link[graphics]{smoothScatter}}?
#' @param loess wartość logiczna - czy do wykresów rozrzutu dodawać linię regresji nieparametrycznej (estymowanej funkcją \code{\link[stats]{loess}})?
#' @param span parametr sterujący wygładzaniem regresji nieparametrycznej (p. \code{\link[stats]{loess}})
#' @param control parametry regresji nieparametrycznej - lista zwrócona przez funkcję \code{\link[stats]{loess.control}}
#' @param paleta funkcja definiująca paletę kolorów do wykorzystania przy rysowaniu wykresów, typowo wynik wywołania funkcji \code{\link[grDevices]{colorRampPalette}}
#' @param pch jeśli smoothScatter=FALSE, znak punktu używanego do rysowania obserwacji (p. \code{\link[graphics]{points}})
#' @param cex jeśli smoothScatter=FALSE, wielkość punktów używanych do rysowania obserwacji (p. \code{\link[graphics]{points}})
#' @param lwd jeśli loess=TRUE, grubość rysowanej linii  (p. \code{\link[graphics]{lines}})
#' @param alpha jeśli smoothScatter=FALSE i nie odpowiada Ci wynik działania domyślnego argumentu ustalajacego, jak bardzo przeźroczyste są punkty, możesz zadać przeźroczystość sam
#' @return funkcja nic nie zwraca
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{\link[lme4]{lmer}}, \code{\link[graphics]{smoothScatter}}, \code{\link[stats]{loess}}
#' @examples
#' x=rnorm(100000)
#' y=x+rnorm(100000)
#' diagnostyka(lm(y ~ x), tytul="y~x")
#' 
#' library(lme4)
#' fm1=lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' diagnostyka(fm1)
#' @export
diagnostyka=function(model, zapiszPng=NULL, tytul=NULL, smoothScatter=TRUE, loess=TRUE, span=0.5, control=loess.control(surface="interpolate", statistics="approximate", trace.hat="approximate", cell=1), paleta=colorRampPalette(c("white", blues9)), pch=16, cex=1, lwd=2, alpha=NULL) {
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
	if ("try-error"%in%c(class(reszty), class(przewidywania))) stop("Na argumencie 'model' musi dać się wykonać funkcje resid() i fitted().")
	
	katalog = getwd()
	if (is.character(zapiszPng)) {
		zapiszPng = try(setwd(zapiszPng[1]), silent=TRUE)
		if ("try-error" %in% class(zapiszPng)) {
			zapiszPng = NULL;
			warning("Nie udało się zapisać wykresów - podano niepoprawną nazwę (ścieżkę) katalogu.")
		}
	}
	parametryGraficzne = par(no.readonly=TRUE)
	par(bg="white")
	
	if (is.null(tytul)) tytul = ""
	podtytul = paste0(paste0(deparse(formula(model), 70), collapse="\n"))
	
	
	if (any(class(model) %in% c("lmerMod", "glmerMod"))) {
		stdReszty = resid(model, scaled=TRUE)
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
			MoreArgs=list(przew=zagrPrzewSt),
			SIMPLIFY=FALSE
		)
		names(zagrPrzewSt) = names(blup)
		efLos = VarCorr(model)  # przydzadzą nam się też odchylenia standardowe efektów losowych
		efLos = lapply(efLos, function(x) return(attributes(x)$stddev))
	} else {
		stdReszty = rstandard(model)
	}
	pierwBezwzglStdReszty = sqrt(abs(stdReszty))
	
	kolorKropek =rgb2hsv(col2rgb(paleta(4)[3]))
	if (is.null(alpha)) alpha = ifelse(
		length(reszty)>100,
		1/log(length(reszty)/10,10),
		1
	)
	kolorKropek = hsv(kolorKropek[1], kolorKropek[2], kolorKropek[3], alpha)
	
	limPrzewidywania         = range(przewidywania)
	limReszty                = max(abs(reszty)) * c(-1,1)
	limStdReszty             = max(abs(stdReszty)) * c(-1,1)
	limPierwBezwzglStdReszty = c(0, max(pierwBezwzglStdReszty))
	# reszty ~ przewidywania
	plot(NA, NA, xlim=limPrzewidywania, ylim=limReszty, main="", xlab="wartości przewidywane", ylab="reszty")
	title(main=tytul, line=2.2)
	title(main=list(podtytul, cex=0.8, font=1), line=0.5)
	if (smoothScatter) {
		smoothScatter(przewidywania, reszty, add=TRUE, nbin=256, bandwidth=c((limPrzewidywania[2] - limPrzewidywania[1]) / 256, (limReszty[2] - limReszty[1]) / 256), colramp=paleta)
	} else {
		points(przewidywania, reszty, pch=pch, cex=cex, col=kolorKropek)
	}
	grid(col=grey(0.5))
	abline(h=0)
	if (loess) {
		lines(przew_npar(reszty ~ przewidywania, data.frame(reszty, przewidywania)), col=2, lwd=lwd, lty=1)
	}
	if (is.character(zapiszPng)) {
		dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_", "reszty_w_funkcji_przewidywania", ".png"), height=1004, width=1004, pointsize=12, res=150)
	}
	# normalność reszt
	temp = qqnorm(stdReszty, plot.it=FALSE)
	plot(NA, NA, xlim=limStdReszty, ylim=limStdReszty, main="", xlab="", ylab="rozkład empiryczny reszt standaryzowanych")
	title(main=tytul, line=2.2)
	title(main=list(podtytul, cex=0.8, font=1), line=0.5)
	title(xlab="rozkład teoretyczny: N(0;1)", line=2)
	title(sub=list("poniżej przekątnej - częstość większa niż teoretyczna,\npowyżej przekątnej - częstość mniejsza niż teoretyczna", cex=0.8, font=1), line=3.8)
	if (smoothScatter) {
		smoothScatter(temp, add=TRUE, transformation=function(x) x^(1/8), nbin=256, bandwidth=c((limStdReszty[2] - limStdReszty[1]) / 512, (limStdReszty[2] - limStdReszty[1]) / 512), colramp=paleta)
	} else {
		points(temp, pch=pch, cex=cex, col=kolorKropek)
	}
	grid(col=grey(0.5))
	abline(0,1, col=2, lwd=2, lty=1)
	if (is.character(zapiszPng)) {
		dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_", "normalnosc_reszt", ".png"), height=1004, width=1004, pointsize=12, res=150)
	}
	# homoscedatyczność
	plot(NA, NA, xlim=limPrzewidywania, ylim=limPierwBezwzglStdReszty, main="", xlab="wartości przewidywane", ylab="")
	title(main=tytul, line=2.2)
	title(main=list(podtytul, cex=0.8, font=1), line=0.5)
	title(ylab=list(as.expression(substitute(sqrt(abs(R)), list(R=as.name("reszty standaryzowane"))))), line=2.5)
	if (smoothScatter) {
		smoothScatter(przewidywania, pierwBezwzglStdReszty, add=TRUE, nbin=256, bandwidth=c((limPrzewidywania[2] - limPrzewidywania[1]) / 256, (limPierwBezwzglStdReszty[2] - limPierwBezwzglStdReszty[1]) / 256), colramp=paleta)
	} else {
		points(przewidywania, pierwBezwzglStdReszty, pch=pch, cex=cex, col=kolorKropek)
	}
	grid(col=grey(0.5))
	if (loess) {
		lines(przew_npar(pierwBezwzglStdReszty~przewidywania, data.frame(pierwBezwzglStdReszty, przewidywania)), col=2, lwd=lwd, lty=1)
	}
	if (is.character(zapiszPng)) {
		dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_", "homoscedatycznosc", ".png"), height=1004, width=1004, pointsize=12, res=150)
	}
	# diagnostyka efektów losowych (z poziomów innych niż indywidualny)
	if (any(class(model) %in% c("lmerMod", "glmerMod"))) {
		for (i in 1:length(blup)) {
			for (j in 1:ncol(blup[[i]])) {
				blupTemp = blup[[i]][[j]]
				zagrPrzewStTemp = zagrPrzewSt[names(blup)[[i]]][[1]]$sr_przew
				# analiza normalności rozkładów BLUPsów;
				skala = unlist(lapply(efLos[grepl(paste0("^", names(blup)[i], "(|[.][[:digit:]]+)$"), names(efLos))], function(x, z) return(x[names(x) == z]), z=names(blup[[i]])[j]))
				temp = qqnorm(blupTemp / skala, plot.it=FALSE)
				limBlup = max(abs(unlist(temp))) * c(-1,1)
				plot(NA, NA, xlim=limBlup, ylim=limBlup, main="", xlab="", ylab=paste0("rozkł. emp. efektu losowego: ", names(blup[[i]])[j], " | ", names(blup)[i]))
				title(main=tytul, line=2.2)
				title(main=list(podtytul, cex=0.8, font=1), line=0.5)
				title(xlab="rozkład teoretyczny: N(0;1)", line=2)
				title(sub=list("poniżej przekątnej - częstość większa niż teoretyczna,\npowyżej przekątnej - częstość mniejsza niż teoretyczna", cex=0.8, font=1), line=3.8)
				if (smoothScatter) {
					smoothScatter(temp, add=TRUE, transformation=function(x) x^(1/8), nbin=256, bandwidth=c((limBlup[2]-limBlup[1]) / 512, (limBlup[2]-limBlup[1]) / 512), colramp=paleta)
				} else {
					points(temp, pch=pch, cex=cex, col=kolorKropek)
				}
				grid(col=grey(0.5))
				abline(0,1, col=2, lwd=2, lty=1)
				if (is.character(zapiszPng)) {
					dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_", "normalnosc_ef_los_", i, "_", j, ".png"), height=1004, width=1004, pointsize=12, res=150)
				}
				# analiza związków między średnimi przewidywanami wynikającymi z części stałej modelu a BLUPsami;
				limBlup = max(abs(unlist(blupTemp))) * c(-1,1)
				limPrzewidywaniaZagr = range(zagrPrzewStTemp)
				if (limPrzewidywaniaZagr[2] > limPrzewidywaniaZagr[1]) {  # dla wyników badań eksperymentalnych ze ślicznie zrównoważonymi schematami może być tak, że nie ma tu czego sprawdzać
					plot(NA, NA, xlim=limPrzewidywaniaZagr, ylim=limBlup, main="", xlab=paste0("śr. wartości przewidywane w ramach grup\nwyróżnionych ze względu na zmienną: ", names(blup)[i]), ylab=paste0("przew. dla efektu losowego: ", names(blup[[i]])[j], " | ", names(blup)[i]))
					title(main=tytul, line=2.2)
					title(main=list(podtytul, cex=0.8, font=1), line=0.5)
					if (smoothScatter) {
						smoothScatter(zagrPrzewStTemp, blupTemp, add=TRUE, nbin=256, bandwidth=c((limPrzewidywaniaZagr[2]-limPrzewidywaniaZagr[1]) / 256, (limBlup[2]-limBlup[1]) / 256), colramp=paleta)
					} else {
						points(zagrPrzewStTemp, blupTemp, pch=pch, cex=cex, col=kolorKropek)
					}
					grid(col=grey(0.5))
					abline(h=0)
					if (loess) {
						lines(przew_npar(blupTemp ~ zagrPrzewStTemp, data.frame(blupTemp, zagrPrzewStTemp)), col=2, lwd=lwd, lty=1)
					}
					if (is.character(zapiszPng)) {
						dev.print(png, paste0(gsub("[\n ]", "_", gsub("[:.]", "", tytul)), "_", "przew_ef_los_w_funkcji_przewidywania", i, "_", j, ".png"), height=1004, width=1004, pointsize=12, res=150)
					}
				} else {
					message(paste0("Dla pogrupowania ze względu na zmienną '", names(blup)[i], "' nie wygenerowano wykresu\n   przewidywania ef. los. | przewidywanie z ef. stałych\nze względu na dosknale zrównoważony schemat analizy (i w efekcie brak zróżnicowania przewidywania z ef. stałych w ramach grup)."))
				}
			}
		}
	}
	# kończenie
	par(parametryGraficzne)
	setwd(katalog)
	invisible(NULL)
}
#' @title Skala kolorow dla wektora liczbowego.
#' @description
#' Funkcja tworzy mapowanie wartości wektora liczbowego na kolory (skalę kolorów) i zwraca wynik takiego mapowania.
#' @param x wektor liczbowy
#' @param kolory wektor tekstowy z nazwami kolorów
#' @param zakres opcjonalny dwualementowy wektor podający zakres wartości, w jakim ma różnicować się skala kolorów
#' @return wektor z nazwami kolorów o długości równej wektorowi x
#' @seealso \code{\link[grDevices]{colorRampPalette}}
#' @examples
#' x=runif(150)
#' plot(x, pch=21, bg=kolorki(x, c("green", "red")))
#' plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red")))
#' plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red"), c(0.4, 0.6)))
#' plot(x, pch=21, bg=kolorki(x, c("green", "yellow", "red"), c(-1, 2)))
#' @export
kolorki=function(x, kolory, zakres=NULL) {
	stopifnot(is.numeric(x), is.character(kolory), is.null(zakres) | is.numeric(zakres))
	stopifnot(length(kolory) > 1, is.null(zakres) | length(zakres) == 2)
	
	if (is.null(zakres)) {  # potnij x na 100 przedziałów o równej szerokości
		x = cut(x, c(-Inf, seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=101)[-c(1,101)], Inf))
	} else {  # potnij x wg przedziałów powstałych przez podział podanego zakresu na 100 kawałków o równej szerokości
		x = cut(x, c(-Inf, seq(zakres[1], zakres[2], length.out=101)[-c(1,101)], Inf))
	}
	paleta = colorRampPalette(kolory)
	lKolorow = length(unique(x))
	return(paleta(nlevels(x))[x])
}
#' @title Konwersja nazw OKE
#' @description
#' Funkcja konwertuje nazwy OKE pomiędzy dwoma standardami stosowanymi w zbiorach EWD - 1) zapisem z użyciem polskich znaków i wielkich liter, 2) zapisem tylko z wykorzystaniem małych liter, bez polskich znaków.
#' @param x czynnik opisujący przydział do OKE
#' @param naPolskie wartość logiczna - czy konwertować do standardu zapisu z polskimi znakami? (jeśli FALSE - konwersja do standardu zapisu bez polskich znaków)
#' @return czynnik
#' @examples
#' x=c("gdansk", "lodz")
#' y=c("Warszawa", "Kraków")
#' nazwy_oke(factor(x))
#' nazwy_oke(factor(y), FALSE)
#' nazwy_oke(factor(c(x, y)))
#' @export
nazwy_oke=function(x, naPolskie=TRUE) {
	stopifnot(is.factor(x), is.logical(naPolskie))
	stopifnot(length(naPolskie) == 1, naPolskie %in% c(TRUE, FALSE))
	stopifnot(all(levels(x) %in% c("gdansk", "jaworzno", "krakow", "lomza", "lodz", "poznan", "warszawa", "wroclaw", "Gdańsk", "Jaworzno", "Kraków", "Łomża", "Łódź", "Poznań", "Warszawa", "Wrocław")))
	if (!(all(levels(x) %in% c("gdansk", "jaworzno", "krakow", "lomza", "lodz", "poznan", "warszawa", "wroclaw")) | all(levels(x) %in% c("Gdańsk", "Jaworzno", "Kraków", "Łomża", "Łódź", "Poznań", "Warszawa", "Wrocław")))) {
		warning("Wykryto poziomy zarówno w formacie z polskimi znakami jak i w formacie bez polskich znaków.")
	}
	if (naPolskie) {
		levels(x)=sub("gdansk",   "Gdańsk",   levels(x))
		levels(x)=sub("jaworzno", "Jaworzno", levels(x))
		levels(x)=sub("krakow",   "Kraków",   levels(x))
		levels(x)=sub("lomza",    "Łomża",    levels(x))
		levels(x)=sub("lodz",     "Łódź",     levels(x))
		levels(x)=sub("poznan",   "Poznań",   levels(x))
		levels(x)=sub("warszawa", "Warszawa", levels(x))
		levels(x)=sub("wroclaw",  "Wrocław",  levels(x))
		x=factor(as.character(levels(x))[x], levels=c("Gdańsk", "Jaworzno", "Kraków", "Łomża", "Łódź", "Poznań", "Warszawa", "Wrocław"))
	} else {
		levels(x)=sub("Gdańsk",   "gdansk",   levels(x))
		levels(x)=sub("Jaworzno", "jaworzno", levels(x))
		levels(x)=sub("Kraków",   "krakow",   levels(x))
		levels(x)=sub("Łomża",    "lomza",    levels(x))
		levels(x)=sub("Łódź",     "lodz",     levels(x))
		levels(x)=sub("Poznań",   "poznan",   levels(x))
		levels(x)=sub("Warszawa", "warszawa", levels(x))
		levels(x)=sub("Wrocław",  "wroclaw",  levels(x))
		x=factor(as.character(levels(x))[x], levels=c("gdansk", "jaworzno", "krakow", "lomza", "lodz", "poznan", "warszawa", "wroclaw"))
	}
	return(x)
}
#' @title Wazone kwantyle
#' @description
#' Funkcja wylicza kwantyle uwzględniając nierówne wagi jednostek obserwacji. Sposób wyliczania kwanrtyli odpowiada podejściu \code{type=2} z funkcji \code{\link[stats]{quantile}}.
#' @param x wektor wartości zmiennej
#' @param wagi wektor wag
#' @param probs wektor prawdopodobieństw definiujących kwantyle, które mają zostać wyliczone
#' @param na.rm wartość logiczna - czy pominąć braki danych?
#' @return wektor z wyliczonymi kwantylami
#' @examples
#' x=1:10
#' w=rep(1:2, each=5)
#' kwantyl_wazony(x, w)
#' quantile(c(1:5, rep(6:10, each=2)), type=2)
#' @export
kwantyl_wazony=function(x, wagi, probs=seq(0, 1, 0.25), na.rm=FALSE) {
	stopifnot(is.numeric(x), is.numeric(wagi), is.numeric(probs), is.logical(na.rm))
	stopifnot(all(wagi >= 0), length(wagi) == length(x), all(probs >= 0), all(probs <= 1), length(na.rm) == 1)
	
	temp = data.frame(x=x, wagi=wagi, wagiOdwr=wagi)
	# obsługa ew. braków danych w x
	if (na.rm) {  # gdzie braki, to zignoruj i po prostu policz na pozostałych
		temp = na.omit(temp)
	} else {
		if (with(temp, any(is.na(x[wagi > 0])))) {  # jeśli jakieś obserwacje z niezerowymi wagami mają braki, zwróć brak danych
			return(setNames(rep(NA, length(probs)), probs))
		} else { # jeśli braki dotyczą tylko obserwacji z wagami równymi zero, można je zignorować
			temp = na.omit(temp)
		}
	}
	# samo wyliczanie kwantyli
	temp = temp[order(temp$x), ]
	temp$wagi    [!is.na(temp$wagi    )] =     cumsum(    temp$wagi    [!is.na(temp$wagi    )] )   / sum(temp$wagi    , na.rm=TRUE)
	temp$wagiOdwr[!is.na(temp$wagiOdwr)] = rev(cumsum(rev(temp$wagiOdwr[!is.na(temp$wagiOdwr)]))) / sum(temp$wagiOdwr, na.rm=TRUE)
	wynik = vector()
	for (i in probs) {
		wynik = c(wynik, mean(temp$x[temp$wagi >= i & temp$wagiOdwr >= (1 - i)]))
	}
	return(setNames(wynik, probs))
}
