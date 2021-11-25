all:
	Rscript runtest.R

manual:
	#sudo apt install r-cran-devtools
	#R -e install.packages("devtools")
	R -e "roxygen2::roxygenise(); devtools::build_manual(path='.')"
	#sudo apt-get install libxml2-dev
	#R -e "install.packages('cli')"
	#R -e "install.packages('knitr')"
	
atp2019:
	R -e "atp2019 = read.csv('data/history2019.csv'); usethis::use_data(atp2019, overwrite = T)"

atpOld:
	R -e "atpOld = read.csv('data/historyOld.csv'); usethis::use_data(atpOld, overwrite = T)"

build:
	R -e "roxygen2::roxygenise('.', roclets=c('rd', 'collate', 'namespace')); devtools::build_manual(path='.')"
	cd ..; R CMD build TrueSkillThroughTime.R
	echo 'importFrom("methods", "new")' >> NAMESPACE
	echo 'importFrom("stats", "dnorm", "pnorm", "qnorm")' >> NAMESPACE
	cd ..; R CMD check --as-cran TrueSkillThroughTime_0.1.0.tar.gz
	
