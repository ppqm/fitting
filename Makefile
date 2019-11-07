
all:
	echo

env:
	echo

fit:
	@pyprofile python fitter/fit.py -p parameters/parameters.pm3.json --method pm3

query:
	@python fitter/fit.py -p parameters/parameters.pm3.test.json2 --method pm3

test:
	@python fitter/mndo.py

