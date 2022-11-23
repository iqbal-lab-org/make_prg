test:
    poetry run pytest

lint:
    poetry run flake8 --max-line-length 88 .

fmt:
    poetry run black .
    poetry run isort .

install:
    poetry install

install-ci:
	python -m pip install --upgrade pip
	python -m pip install poetry
	poetry install --no-interaction

check-fmt:
    poetry run black --check .
    poetry run isort --check-only .

precommit: fmt lint test

build:
    poetry run build