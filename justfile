PROJECT := "make_prg"
OPEN := if os() == "macos" { "open" } else { "xdg-open" }
VERSION := `poetry version -s`

# run tests (locally)
test *FLAGS:
    poetry run pytest {{FLAGS}}

# run tests with coverage and open in a browser
coverage:
    poetry run pytest --cov-report html
    {{ OPEN }} htmlcov/index.html

# lint code with flake8
lint:
    poetry run flake8 . --extend-exclude=venv,.venv

# format code with black and isort
fmt:
    poetry run black .
    poetry run isort .

install:
    poetry install

install-ci:
    poetry install --no-interaction

# check code formatting with black and isort
check-fmt:
    poetry run black --check .
    poetry run isort --check-only .

# recipes to run before commiting: format, lint, and test
precommit: fmt lint test

# build a python release
build:
    poetry run build

# prints out the commands to run to tag the release and push it
tag:
    @echo "Run \`git tag -a {{ VERSION }} -m <message>\` to tag the release"
    @echo "Then run \`git push origin {{ VERSION }}\` to push the tag"
