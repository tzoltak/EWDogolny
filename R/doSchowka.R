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
				daneM = data.frame(x$model[,1], model.matrix(x)[, -1], check.names=FALSE)
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
				odchStd = unlist(lapply(as.data.frame(model.matrix(x)[, -1]), sd))
				# standaryzacja idzie względem wariancji zm. zależnej na poziomie indywidualnym (wariancja przewidywania + wariancja błędów indywidualnych)
				warPrzew = var(predict(x, re.form=~0), na.rm=TRUE)
				odchStd = c(zal=sqrt(warPrzew + sigma(x)^2), odchStd)
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
