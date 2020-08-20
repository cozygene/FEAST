.PHONY: test

test:
	nosetests -v -s q2_feast --with-coverage --cover-package=q2_feast
