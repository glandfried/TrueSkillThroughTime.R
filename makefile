all:
	Rscript runtest.R

manual:
	#sudo apt install r-cran-devtools
	R -e "roxygen2::roxygenise(); devtools::build_manual(path='.')"
	
atp2019:
	R -e "atp2019 = read.csv('data/history2019.csv'); usethis::use_data(atp2019, overwrite = T)"

atpOld:
	R -e "atpOld = read.csv('data/historyOld.csv'); usethis::use_data(atpOld, overwrite = T)"

build: manual
	cd ..; R CMD build TrueSkillThroughTime.R
	cd ..; R CMD check --as-cran TrueSkillThroughTime_0.1.0.tar.gz
	
