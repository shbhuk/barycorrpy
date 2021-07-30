init: 
	python setup.py develop
	pip install pytest pytest-coverage coveralls
docker:
	docker build --tag bary:latest .
	docker run -it -v ${PWD}:/code/barycorrpy bary:latest bash

regression_tests:
	pytest --cov=barycorrpy --pyargs barycorrpy.tests
	coveralls
