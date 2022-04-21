/*

Nombre:      ardl
Descripción: Do file para el documento de investigación "Pronósticos Uniecuacionales 
			 para los Ingresos Tributarios de Honduras". Publicado en la Revista de 
			 Administración Tributaria, N° 46, octubre 2020.
			 https://biblioteca.ciat.org/opac/book/5743
Autor:       Jose Carlo Bermúdez, jbermudez@sar.gob.hn

Departamento de Estudios Fiscales y Económicos
Direccion Nacional de Gestion Estratégica
Servicio de Administracion de Rentas, Honduras

*/

	clear all 
	clear matrix 
	global path "C:\Users\jbermudez\OneDrive - SAR\Notas técnicas y papers\Modelo ARDL" 	// Ingrese directorio 
	cd "$path"
	use data_ardl.dta, clear

	tsmktim tiempo, start(2007m1) //Definicion de variable de tiempo y dummies estacionales (Definition of time series variable and seasonal dummies)
	tsset tiempo
	generate m=month(dofm(tiempo))
	tabulate m, generate(dum)
	gen trend=_n

	foreach var of varlist imae ener itrib imae_sa ener_sa itrib_sa {  //Rescalamiento de variables (Rescaling variables)
		gen ln_`var'=log(`var')
	}


	* Gráfico 3 de las variables
	foreach var of varlist ln_imae ln_itrib ln_ener {
		tsline `var' `var'_sa if tiempo<=tm(2019m12), ytitle("") xtitle("") graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(blue red)
		graph export "$path\graph3_`var'.pdf", replace
	}


	* Gráfico 1. Rolling Correlations
	set seed 2000
	
	mvcorr ln_itrib_sa ln_imae_sa, win(12) gen(rho1) 
	mvcorr ln_itrib_sa ln_ener_sa, win(12) gen(rho2)

	bootstrap se_rho1 = (r(sd)/sqrt(145)), reps(400): summ rho1 		//Error estándar mediante bootstrapping 
	gen rho1_upper= rho1 + (1.96 * e(b_bs)[1,1])				 				
	gen rho1_lower= rho1 - (1.96 * e(b_bs)[1,1])
	
	bootstrap se_rho2 = (r(sd)/sqrt(145)), reps(400): summ rho2
	gen rho2_upper= rho2 + (1.96 * e(b_bs)[1,1])
	gen rho2_lower= rho2 - (1.96 * e(b_bs)[1,1])

	graph twoway rarea rho1_upper rho1_lower tiempo if tiempo>=tm(2008m1), color(blue*.5) || ///
	             line rho1 tiempo if tiempo>=tm(2008m1), lpattern(dash) lwidth(medthick) lcolor(blue*1.25)||, ///
				 ylabel(-1(0.5)1) ytitle("") xtitle("") graphregion(fcolor(white)) ///
				 legend(row(2) ring(0) position(4) label(1 "± 95%") label(2 "Rho 1") order(2 1))
				 graph export "$path\graph_1a.pdf", replace
				 
	graph twoway rarea rho2_upper rho2_lower tiempo if tiempo>=tm(2008m1), color(red*.5) || ///
	             line rho2 tiempo if tiempo>=tm(2008m1), lpattern(dash) lwidth(medthick) lcolor(red*1.25)||, ///
				 ylabel(-1(0.5)1) ytitle("") xtitle("") graphregion(fcolor(white)) ///
				 legend(row(2) ring(0) position(4) label(1 "± 95%") label(2 "Rho 2") order(2 1))
				 graph export "$path\graph_1b.pdf", replace


	* Prueba de raíz unitaria empleando Dolado et al (1990)
	varsoc ln_imae_sa, exog(tiempo)                
	dfuller ln_imae_sa, lags(3) trend reg
	dfuller ln_imae_sa, lags(3) reg
	dfuller ln_imae_sa, lags(3) noconstant reg

	varsoc ln_itrib_sa, exog(tiempo) 				
	dfuller ln_itrib_sa, lags(2) trend reg
	dfuller ln_itrib_sa, lags(2) reg
	dfuller ln_itrib_sa, lags(2) noconstant reg

	varsoc ln_ener_sa, exog(tiempo)					
	dfuller ln_ener_sa, lags(1) trend reg
	dfuller ln_ener_sa, lags(1) reg
	dfuller ln_ener_sa, lags(1) noconstant reg 


	*	Modelo ARDL mediante metodología de cointegración a la Pesaran, Shin & Smith (2001)
	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) aic dots  							//Estimacion de rezagos óptimos 
	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) ec1 lags(3 3 2) regstore(ardl_ec) 	//ARDL (p,q,r) en forma de MCE (Tabla 2 - ecuación 2)
	estat ectest 																									//Verifico Cointegración mediante prueba de Pesaran et al (2001)
	estimates restore ardl_ec
	`e(cmdline)' vce(robust)
										
	predict ardl_1, residuals 						//Reservo los residuos del MCE (saving ECM residuals)
	estimates restore ardl_ec
	estat bgodfrey, lags(1) small					//Autocorrelación (checking for autocorrelation)
	estat hettest									//Homocedasticidad (checking for heterokedasticity)
	estat imtest, white

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl)			// Especificación Modelo ARDL Óptimo (Tabla 2 - ecuación 3)
	estimates restore ardl
	`e(cmdline)' vce(robust)

	* ------------------------------------------------------------------
	* Estimación de diferentes especificaciones de modelos univariados
	* ------------------------------------------------------------------
	forval p=1/12{
		qui arima ln_itrib, ar(`p') nolog //Modelo AR(p) óptimo
			estimates store ar_`p'
	}
	estimates stats ar_* 			 //Sugiere un AR(3)

	forval q=1/12{
		qui arima ln_itrib, ma(`q') nolog //Modelo MA(q) óptimo
			estimates store ma_`q'
	}
	estimates stats ma_* 			 //Sugiere un MA(8)

	forval p=1/4{
		forval q=1/4{
			qui arima ln_itrib, ar(`p') ma(`q') nolog //Modelo ARMA(p,q) óptimo
				estimates store arma_`p'_`q'
		}
	}
	estimates stats arma_*_*   			//Sugiere un ARMA(3,3)

	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib ln_imae if tin(, 2019m12), ar(1/`p') ma(1/`q') nolog //Modelo ARMAX(p,q) con el IMAE
				estimates store armaxi_`p'_`q'
		}
	}
	estimates stats armaxi_*_*			//Sugiere un ARMAX(2,1)

	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib ln_ener if tin(, 2019m12), ar(1/`p') ma(1/`q') nolog //Modelo ARMAX(p,q) con la Energía
				estimates store armaxe_`p'_`q'
		}
	}
	estimates stats armaxe_*_*			//Sugiere un ARMAX(1,1)

	set matsize 2000
	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib, arima(`p',1,0) sarima(`p',1,0,12) nolog //Modelo SAR(p,d,0)(P,D,Q)
				estimates store sar_`p'_1_0_`p'_1_0
			}
	}
	estimates stats sar_*_1_0_*_1_0		//Sugiere un SAR(3,0)(3,0)
		
	set matsize 2000
	forval p=1/3{
		forval q=1/3{
			qui arima ln_itrib, arima(`p',1,`q') sarima(`p',1,`q',12) nolog //Modelo SARIMA(p,d,q)(P,D,Q)
				estimates store sarima_`p'_1_`q'_`p'_1_`q'
			}
		}
	estimates stats sarima_*_1_*_*_1_*	//Sugiere un SARIMA(2,1)(2,1)
	

	* ------------------------------------------------------------------
	* 	Análisis de pronóstico dentro de muestra h=1
	* ------------------------------------------------------------------
	qui arima ln_itrib, ar(3) vce(robust) 
	predict fitted_ar 

	qui arima ln_itrib, ma(8) vce(robust)
	predict fitted_ma

	qui arima ln_itrib, ar(3) ma(3) vce(robust)
	predict fitted_arma

	qui arima ln_itrib ln_imae, ar(2) ma(1) vce(robust)
	predict fitted_armaxi

	qui arima ln_itrib ln_ener, ar(1) ma(1) vce(robust)
	predict fitted_armaxe

	forval p=1/3 {
		forval q=1/2 {
			quietly reg ln_itrib L(1/`p').ln_itrib L(1/`p').ln_ener L(1/`q').ln_imae dum1 dum4 dum6 dum8 dum10 dum12 trend, vce(robust) 
				predict ardl_fitted_`p'_`p'_`q'
		}
	}

	matrix drop _all
	foreach var of varlist fitted_* ardl_fitted_*_* {
		fcstats ln_itrib `var'   								//Análisis de precisión pronósticos para h=1
		mat forecast_h1 = nullmat(forecast_h1) \ (r(rmse), r(mae), r(mape), r(theilu))
	}
	
	* Tabla 3
	frmttable using "$path/fcst_h1", replace tex statmat(forecast_h1) ctitle("Modelo", "RECM", "EAM", "EPAM", "Theil") fr basefont(tiny) sd(4,4,4,4) ///
			  rtitles("AR(3)"\"MA(8)"\"ARMA(3,3)"\"ARMAX (2,1)"\"ARMAX (1,1)"\"ARDL (1,1,1)"\"ARDL (1,1,2)"\"ARDL (2,2,1)"\"ARDL (2,2,2)"\"ARDL (3,3,1)"\"ARDL (3,3,2)")	          
	
	* Gráfico 4
	local fitted "ardl_fitted_3_3_2 fitted_ar fitted_ma fitted_arma fitted_armaxi fitted_armaxe ardl_fitted_1_1_1 ardl_fitted_1_1_2 ardl_fitted_2_2_1 ardl_fitted_2_2_2 ardl_fitted_3_3_1"
	tsline ln_itrib `fitted' if tiempo>=tm(2015m1), ylabel(8.2(0.5)9.8) ytitle("") xtitle("") ///
	       graphregion(fcolor(white)) lwidth(medthick medthick thick thick thick thick thick thick thick thick thick thick) ///
		   lcolor(blue red ltblue ltblue ltblue ltblue ltblue ltblue ltblue ltblue ltblue ltblue) ///
		   lpattern(solid dash dot dot dot dot dot dot dot dot dot dot) ///
		   legend(order(1 "Observado" 2 "ARDL(3,3,2)" 3 "Ajustados"))
	graph export "$path\graph4.pdf", replace
	
	drop fitted_* ardl_fitted_*_*


	* ------------------------------------------------------------------
	* Análisis de pronóstico fuera de muestra h=3 
	* ------------------------------------------------------------------
	qui arima ln_itrib if tiempo<=tm(2019m9), ar(3) vce(robust) 								
	estimates store ar_3h
	forecast create ar_model_3h, replace
	forecast estimates ar_3h
	forecast solve, suffix(_ar_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m9), ma(8) vce(robust)									
	estimates store ma_3h
	forecast create ma_model_3h, replace
	forecast estimates ma_3h
	forecast solve, suffix(_ma_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m9), ar(3) ma(3) vce(robust) 							
	estimates store arma_3h
	forecast create arma_model_3h, replace
	forecast estimates arma_3h
	forecast solve, suffix(_arma_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib ln_imae if tiempo<=tm(2019m9), ar(2) ma(1) vce(robust)					
	estimates store armaxi_3h
	forecast create armaxi_model_3h, replace
	forecast estimates armaxi_3h
	forecast solve, suffix(_armaxi_3h) begin(tm(2019m10)) log(off)

	qui arima ln_itrib ln_ener if tiempo<=tm(2019m9), ar(1) ma(1) vce(robust)      				
	estimates store armaxe_3h
	forecast create armaxe_model_3h, replace
	forecast estimates armaxe_3h
	forecast solve, suffix(_armaxe_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl_1_1_1) 
	estimates restore ardl_1_1_1
	`e(cmdline)' vce(robust)
	estimates store ardl_1
	forecast create ardl_1_model_3h, replace
	forecast estimates ardl_1
	forecast solve, suffix(_ardl_1_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl_1_1_2) 
	estimates restore ardl_1_1_2
	`e(cmdline)' vce(robust)
	estimates store ardl_2
	forecast create ardl_2_model_3h, replace
	forecast estimates ardl_2
	forecast solve, suffix(_ardl_2_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl_2_2_1) 
	estimates restore ardl_2_2_1
	`e(cmdline)' vce(robust)
	estimates store ardl_3
	forecast create ardl_3_model_3h, replace
	forecast estimates ardl_3
	forecast solve, suffix(_ardl_3_3h) begin(tm(2019m10)) log(off)

	quietly ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl_2_2_2) 
	estimates restore ardl_2_2_2
	`e(cmdline)' vce(robust)
	estimates store ardl_4
	forecast create ardl_4_model_3h, replace
	forecast estimates ardl_4
	forecast solve, suffix(_ardl_4_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl_3_3_1) 
	estimates restore ardl_3_3_1
	`e(cmdline)' vce(robust)
	estimates store ardl_5
	forecast create ardl_5_model_3h, replace
	forecast estimates ardl_5
	forecast solve, suffix(_ardl_5_3h) begin(tm(2019m10)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m9), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl_3_3_2)
	estimates restore ardl_3_3_2
	`e(cmdline)' vce(robust)
	estimates store ardl_6
	forecast create ardl_6_model_3h, replace
	forecast estimates ardl_6
	forecast solve, suffix(_ardl_6_3h) begin(tm(2019m10)) log(off)

	matrix drop _all
	foreach var of varlist ln_itrib_ar_3h ln_itrib_ma_3h ln_itrib_arma_3h ln_itrib_armaxi_3h ///
	                       ln_itrib_armaxe_3h ln_itrib_ardl_1_3h ln_itrib_ardl_2_3h ///
						   ln_itrib_ardl_3_3h ln_itrib_ardl_4_3h ln_itrib_ardl_5_3h ln_itrib_ardl_6_3h {
		fcstats ln_itrib `var'   													//Análisis de precisión pronósticos para h=3	
		mat forecast_h3 = nullmat(forecast_h3) \ (r(rmse), r(mae), r(mape), r(theilu))
	}
	
		* Tabla 4
	frmttable using "$path/fcst_h3", replace tex statmat(forecast_h3) ctitle("Modelo", "RECM", "EAM", "EPAM", "Theil") fr basefont(tiny) sd(4,4,4,4) ///
			  rtitles("AR(3)"\"MA(8)"\"ARMA(3,3)"\"ARIMAX (2,1,1)"\"ARIMAX (1,1,1)"\"ARDL (1,1,1)"\"ARDL (1,1,2)"\"ARDL (2,2,1)"\"ARDL (2,2,2)"\"ARDL (3,3,1)"\"ARDL (3,3,2)")	          

	drop ln_itrib_*_* ln_itrib_*_*_*

	* ------------------------------------------------------------------
	* Análisis de pronóstico fuera de muestra h=6 
	* ------------------------------------------------------------------
	qui arima ln_itrib if tiempo<=tm(2019m6), ar(3) vce(robust) 								
	estimates store ar_6h
	forecast create ar_model_6h, replace
	forecast estimates ar_6h
	forecast solve, suffix(_ar_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m6), ma(8) vce(robust)									
	estimates store ma_6h
	forecast create ma_model_6h, replace
	forecast estimates ma_6h
	forecast solve, suffix(_ma_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib if tiempo<=tm(2019m6), ar(3) ma(3) vce(robust) 							
	estimates store arma_6h
	forecast create arma_model_6h, replace
	forecast estimates arma_6h
	forecast solve, suffix(_arma_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib ln_imae if tiempo<=tm(2019m6), ar(2) ma(1) vce(robust)					
	estimates store armaxi_6h
	forecast create armaxi_model_6h, replace
	forecast estimates armaxi_6h
	forecast solve, suffix(_armaxi_6h) begin(tm(2019m7)) log(off)

	qui arima ln_itrib ln_ener if tiempo<=tm(2019m6), ar(1) ma(1) vce(robust)      				
	estimates store armaxe_6h
	forecast create armaxe_model_6h, replace
	forecast estimates armaxe_6h
	forecast solve, suffix(_armaxe_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl1_1_1_1) 
	estimates restore ardl1_1_1_1
	`e(cmdline)' vce(robust)
	estimates store ardl1_1
	forecast create ardl_1_model_6h, replace
	forecast estimates ardl1_1
	forecast solve, suffix(_ardl_1_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl1_1_1_2)
	estimates restore ardl1_1_1_2
	`e(cmdline)' vce(robust)
	estimates store ardl1_2
	forecast create ardl_2_model_6h, replace
	forecast estimates ardl1_2
	forecast solve, suffix(_ardl_2_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl1_2_2_1) 
	estimates restore ardl1_2_2_1
	`e(cmdline)' vce(robust)
	estimates store ardl1_3
	forecast create ardl_3_model_6h, replace
	forecast estimates ardl1_3
	forecast solve, suffix(_ardl_3_6h) begin(tm(2019m7)) log(off)

	quietly ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl1_2_2_2) 
	estimates restore ardl1_2_2_2
	`e(cmdline)' vce(robust)
	estimates store ardl1_4
	forecast create ardl_4_model_6h, replace
	forecast estimates ardl1_4
	forecast solve, suffix(_ardl_4_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl1_3_3_1) 
	estimates restore ardl1_3_3_1
	`e(cmdline)' vce(robust)
	estimates store ardl1_5
	forecast create ardl_5_model_6h, replace
	forecast estimates ardl1_5
	forecast solve, suffix(_ardl_5_6h) begin(tm(2019m7)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2019m6), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl1_3_3_2) 
	estimates restore ardl1_3_3_2
	`e(cmdline)' vce(robust)
	estimates store ardl1_6
	forecast create ardl_6_model_6h, replace
	forecast estimates ardl1_6
	forecast solve, suffix(_ardl_6_6h) begin(tm(2019m7)) log(off)

	matrix drop _all
	foreach var of varlist ln_itrib_ar_6h ln_itrib_ma_6h ln_itrib_arma_6h ln_itrib_armaxi_6h ///
					       ln_itrib_armaxe_6h ln_itrib_ardl_1_6h ln_itrib_ardl_2_6h ///
						   ln_itrib_ardl_3_6h ln_itrib_ardl_4_6h ln_itrib_ardl_5_6h ln_itrib_ardl_6_6h {
		fcstats ln_itrib `var'   													//Análisis de precisión pronósticos para h=6
		mat forecast_h6 = nullmat(forecast_h6) \ (r(rmse), r(mae), r(mape), r(theilu))
	}
	
	* Tabla 5
	frmttable using "$path/fcst_h6", replace tex statmat(forecast_h6) ctitle("Modelo", "RECM", "EAM", "EPAM", "Theil") fr basefont(tiny) sd(4,4,4,4) ///
			  rtitles("AR(3)"\"MA(8)"\"ARMA(3,3)"\"ARIMAX (2,1,1)"\"ARIMAX (1,1,1)"\"ARDL (1,1,1)"\"ARDL (1,1,2)"\"ARDL (2,2,1)"\"ARDL (2,2,2)"\"ARDL (3,3,1)"\"ARDL (3,3,2)")	          

	drop ln_itrib_*_* ln_itrib_*_*_*

	* ------------------------------------------------------------------
	* Análisis de pronóstico fuera de muestra h=12 
	* ------------------------------------------------------------------
	qui arima ln_itrib if tiempo<=tm(2018m12), ar(3) vce(robust) 								
	estimates store ar_12h
	forecast create ar_model_12h, replace
	forecast estimates ar_12h
	forecast solve, suffix(_ar_12h) begin(tm(2019m1)) log(off)

	qui arima ln_itrib if tiempo<=tm(2018m12), ma(8) vce(robust)									
	estimates store ma_12h
	forecast create ma_model_12h, replace
	forecast estimates ma_12h
	forecast solve, suffix(_ma_12h) begin(tm(2019m1)) log(off)

	qui arima ln_itrib if tiempo<=tm(2018m12), ar(3) ma(3) vce(robust) 							
	estimates store arma_12h
	forecast create arma_model_12h, replace
	forecast estimates arma_12h
	forecast solve, suffix(_arma_12h) begin(tm(2019m1)) log(off)

	qui arima ln_itrib ln_imae if tiempo<=tm(2018m12), ar(2) ma(1) vce(robust)		
	estimates store armaxi_12h
	forecast create armaxi_model_12h, replace
	forecast estimates armaxi_12h
	forecast solve, suffix(_armaxi_12h) begin(tm(2019m1)) log(off)

	*qui arima ln_itrib ln_ener if tiempo<=tm(2018m12), ar(1) ma(1) vce(robust) *// <- No se estima bien la log likelihood   
	*estimates store armaxe_12h
	*forecast create armaxe_model_12h, replace
	*forecast estimates armaxe_12h
	*forecast solve, suffix(_armaxe_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl2_1_1_1)
	estimates restore ardl2_1_1_1
	`e(cmdline)' vce(robust)
	estimates store ardl2_1
	forecast create ardl_1_model_12h, replace
	forecast estimates ardl2_1
	forecast solve, suffix(_ardl_1_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl2_1_1_2) 
	estimates restore ardl2_1_1_2
	`e(cmdline)' vce(robust)
	estimates store ardl2_2
	forecast create ardl_2_model_12h, replace
	forecast estimates ardl2_2
	forecast solve, suffix(_ardl_2_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl2_2_2_1) 
	estimates restore ardl2_2_2_1
	`e(cmdline)' vce(robust)
	estimates store ardl2_3
	forecast create ardl_3_model_12h, replace
	forecast estimates ardl2_3
	forecast solve, suffix(_ardl_3_12h) begin(tm(2019m1)) log(off)

	quietly ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl2_2_2_2)
	estimates restore ardl2_2_2_2
	`e(cmdline)' vce(robust)
	estimates store ardl2_4
	forecast create ardl_4_model_12h, replace
	forecast estimates ardl2_4
	forecast solve, suffix(_ardl_4_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl2_3_3_1) 
	estimates restore ardl2_3_3_1
	`e(cmdline)' vce(robust)
	estimates store ardl2_5
	forecast create ardl_5_model_12h, replace
	forecast estimates ardl2_5
	forecast solve, suffix(_ardl_5_12h) begin(tm(2019m1)) log(off)

	ardl ln_itrib ln_ener ln_imae if tiempo<=tm(2018m12), exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl2_3_3_2) 
	estimates restore ardl2_3_3_2
	`e(cmdline)' vce(robust)
	estimates store ardl2_6
	forecast create ardl_6_model_12h, replace
	forecast estimates ardl2_6
	forecast solve, suffix(_ardl_6_12h) begin(tm(2019m1)) log(off)

	matrix drop _all
	foreach var of varlist ln_itrib_ar_12h ln_itrib_ma_12h ln_itrib_arma_12h ln_itrib_armaxi_12h ///
						   ln_itrib_ardl_1_12h ln_itrib_ardl_2_12h ln_itrib_ardl_3_12h ///
						   ln_itrib_ardl_4_12h ln_itrib_ardl_5_12h ln_itrib_ardl_6_12h {
		fcstats ln_itrib `var'   													//Análisis de precisión pronósticos para h=12
		mat forecast_h12 = nullmat(forecast_h12) \ (r(rmse), r(mae), r(mape), r(theilu))
	}
	
	* Tabla 6
	frmttable using "$path/fcst_h12", replace tex statmat(forecast_h12) ctitle("Modelo", "RECM", "EAM", "EPAM", "Theil") fr basefont(tiny) sd(4,4,4,4) ///
			  rtitles("AR(3)"\"MA(8)"\"ARMA(3,3)"\"ARIMAX (2,1,1)"\"ARDL (1,1,1)"\"ARDL (1,1,2)"\"ARDL (2,2,1)"\"ARDL (2,2,2)"\"ARDL (3,3,1)"\"ARDL (3,3,2)")	          

	drop ln_itrib_*_* ln_itrib_*_*_* 
	

	* ------------------------------------------------------------------
	* Estimación de multiplicadores de variables endógenas
	* ------------------------------------------------------------------
	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 1) regstore(ardl1) 
	estimates restore ardl1
	`e(cmdline)' vce(robust)
	gen mult_imae1=_b[ln_imae] in 1 												
	replace mult_imae1=L.mult_imae1*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
		forval c=3/12{
			replace mult_imae1=L.mult_imae1*_b[L1.ln_itrib] in `c'
	}
	gen mult_ener1=_b[ln_ener] in 1
	replace mult_ener1=L.mult_ener1*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				
		forval c=3/12{
			replace mult_ener1=L.mult_ener1*_b[L1.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(1 1 2) regstore(ardl2) 
	estimates restore ardl2
	`e(cmdline)' vce(robust)
	gen mult_imae2=_b[ln_imae] in 1 												
	replace mult_imae2=L.mult_imae2*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae2=L.mult_imae2*_b[L1.ln_itrib]+_b[L2.ln_imae] in 3
		forval c=4/12{
			replace mult_imae2=L.mult_imae2*_b[L1.ln_itrib] in `c'
	}
	gen mult_ener2=_b[ln_ener] in 1
	replace mult_ener2=L.mult_ener2*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				
		forval c=3/12{
			replace mult_ener2=L.mult_ener2*_b[L1.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 1) regstore(ardl3) 
	estimates restore ardl3
	`e(cmdline)' vce(robust)
	gen mult_imae3=_b[ln_imae] in 1 												
	replace mult_imae3=L.mult_imae3*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
		forval c=3/12{
			replace mult_imae3=L.mult_imae3*_b[L2.ln_itrib] in `c'
	}
	gen mult_ener3=_b[ln_ener] in 1
	replace mult_ener3=L.mult_ener3*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				
	replace mult_ener3=L.mult_ener3*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
		forval c=4/12{
			replace mult_ener3=L.mult_ener3*_b[L2.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(2 2 2) regstore(ardl4) 
	estimates restore ardl4
	`e(cmdline)' vce(robust)
	gen mult_imae4=_b[ln_imae] in 1 												
	replace mult_imae4=L.mult_imae4*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae4=L.mult_imae4*_b[L2.ln_itrib]+_b[L2.ln_imae] in 3
		forval c=4/12{
			replace mult_imae4=L.mult_imae4*_b[L2.ln_itrib] in `c'
	}
	gen mult_ener4=_b[ln_ener] in 1
	replace mult_ener4=L.mult_ener4*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				
	replace mult_ener4=L.mult_ener4*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
		forval c=4/12{
			replace mult_ener4=L.mult_ener4*_b[L2.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 1) regstore(ardl5) 
	estimates restore ardl5
	`e(cmdline)' vce(robust)
	gen mult_imae5=_b[ln_imae] in 1 												
	replace mult_imae5=L.mult_imae5*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae5=L.mult_imae5*_b[L2.ln_itrib] in 3
		forval c=4/12{
			replace mult_imae5=L.mult_imae5*_b[L3.ln_itrib] in `c'
	}
	gen mult_ener5=_b[ln_ener] in 1
	replace mult_ener5=L.mult_ener5*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				
	replace mult_ener5=L.mult_ener5*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
	replace mult_ener5=L.mult_ener5*_b[L3.ln_itrib]+_b[L3.ln_ener] in 4
		forval c=5/12{
			replace mult_ener5=L.mult_ener5*_b[L3.ln_itrib] in `c'
	}

	ardl ln_itrib ln_ener ln_imae, exog(dum1 dum4 dum6 dum8 dum10 dum12 trend) lags(3 3 2) regstore(ardl6) 
	estimates restore ardl6
	`e(cmdline)' vce(robust)
	gen mult_imae6=_b[ln_imae] in 1 												
	replace mult_imae6=L.mult_imae6*_b[L1.ln_itrib]+_b[L1.ln_imae] in 2
	replace mult_imae6=L.mult_imae6*_b[L2.ln_itrib]+_b[L2.ln_imae] in 3
		forval c=4/12{
			replace mult_imae6=L.mult_imae6*_b[L3.ln_itrib] in `c'
	}
	gen mult_ener6=_b[ln_ener] in 1
	replace mult_ener6=L.mult_ener6*_b[L1.ln_itrib]+_b[L1.ln_ener] in 2				
	replace mult_ener6=L.mult_ener6*_b[L2.ln_itrib]+_b[L2.ln_ener] in 3
	replace mult_ener6=L.mult_ener6*_b[L3.ln_itrib]+_b[L3.ln_ener] in 4
		forval c=5/12{
			replace mult_ener6=L.mult_ener6*_b[L3.ln_itrib] in `c'
	}

	* Gráfico 2a
	local mult_imae "mult_imae1 mult_imae2 mult_imae3 mult_imae4 mult_imae5 mult_imae6"  
	egen mult_imae_mean=rowmean(`mult_imae')
	local y=0
	line mult_imae_mean mult_imae6 trend if trend<=12, ylabel(-.7(0.3).8) ytitle("") yline(`y') xlabel(1(1)12) xtitle("") ///
	     graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(red blue) lpattern(longdash solid) legend(order(1 "Promedio" 2 "ARDL(3,3,2)"))
		 graph export "$path\graph_2a.pdf", replace

	* Gráfico 2b
	local mult_ener "mult_ener1 mult_ener2 mult_ener3 mult_ener4 mult_ener5 mult_ener6"  
	egen mult_ener_mean=rowmean(`mult_ener')
	local y=0
	line mult_ener_mean mult_ener6 trend if trend<=12, ylabel(-.01(0.2).67) ytitle("") yline(`y') xlabel(1(1)12) xtitle("") ///
	     graphregion(fcolor(white)) lwidth(medthick medthick) lcolor(red blue) lpattern(longdash solid) legend(order(1 "Promedio" 2 "ARDL(3,3,2)"))
		 graph export "$path\graph_2b.pdf", replace

