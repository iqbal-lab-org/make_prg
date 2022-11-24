test *FLAGS:
    poetry run pytest {{FLAGS}}

lint:
    poetry run flake8 --max-line-length 88 .

fmt:
    poetry run black .
    poetry run isort .

install:
    poetry install

install-ci:
    poetry config experimental.new-installer false
    poetry install --no-interaction
    poetry env info

check-fmt:
    poetry env info
    poetry run black --check .
    poetry run isort --check-only .

precommit: fmt lint test

build:
    poetry run build