.PHONY: test

test:
	nosetests -v -s q2_FEAST --with-coverage --cover-package=q2_FEAST
