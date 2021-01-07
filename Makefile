PYTHON ?= python
PYTEST ?= pytest
CTAGS ?= ctags

# installation instructions are not clear to me:
#   setup.py install -> the package cannot be imported from outside the repository even though he is 
#   visible in the conda list as <develop> in the same way as when 
#   setup.py develop -> everything is ok
#install: clean-ctags
#	$(PYTHON) setup.py install

install:
	$(PYTHON) setup.py develop 

test:
	$(PYTEST) --cov-config=.coveragerc --cov-report term-missing --cov geneids geneids

ctags:
	$(CTAGS) --python-kinds=-i --exclude=*/tests/* -R geneids

clean:
	rm -f tags
