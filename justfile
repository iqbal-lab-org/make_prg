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
    poetry install --no-interaction

install-mafft:
    if [ {{os()}} = "macos" ]; then wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-signed.pkg -O mafft.pkg && sudo installer -pkg mafft.pkg -target / && mafft --version ; elif [ {{os()}} = "linux" ]; then wget https://mafft.cbrc.jp/alignment/software/mafft_7.505-1_amd64.deb -O mafft.deb && sudo dpkg -i mafft.deb && mafft --version; else echo "Only support installing mafft for linux or mac"; exit 1; fi

check-fmt:
    poetry run black --check .
    poetry run isort --check-only .

precommit: fmt lint test

build:
    poetry run build