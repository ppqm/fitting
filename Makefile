
all:
	echo

env:
	echo

fit:
	@anaconda fitter/fit.py -p parameters/parameters.pm3.json --method pm3

query:
	@anaconda fitter/fit.py -p parameters/parameters.pm3.test.json2 --method pm3

test:
	@anaconda fitter/mndo.py

