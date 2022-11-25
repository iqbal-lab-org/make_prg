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
    poetry run flake8 . --exclude venv

# format code with black and isort
fmt:
    poetry run black . --exclude venv
    poetry run isort . --skip venv

install:
    poetry install

install-ci: install-mafft
    poetry install --no-interaction

# do a local install of MAFFT (intended for use on CI)
install-mafft:
    if [ {{os()}} = "macos" ]; then wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-signed.pkg -O mafft.pkg && sudo installer -pkg mafft.pkg -target / && rm mafft.pkg && mafft --version ; elif [ {{os()}} = "linux" ]; then wget https://mafft.cbrc.jp/alignment/software/mafft_7.505-1_amd64.deb -O mafft.deb && sudo dpkg -i mafft.deb && rm mafft.deb && mafft --version; else echo "Only support installing mafft for linux or mac"; exit 1; fi

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
