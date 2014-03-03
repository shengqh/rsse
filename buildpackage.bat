E:
cd "E:/sqh/programs/rsse"

R --vanilla -f "buildpackage.R"

rm source/inst/doc/rsse-concordance.tex
rm source/inst/doc/rsse.log
rm source/inst/doc/rsse.synctex.gz
rm source/inst/doc/rsse.tex
rm source/inst/doc/rsse.toc
rm source/inst/doc/rsse.bbl

R CMD check source
R CMD build source
R CMD INSTALL --build source
Rscript -e "install.packages('rsse_0.98.2.zip', repos = NULL)"

