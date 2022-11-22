FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y git python3 python3-dev python3-pip graphviz libgraphviz-dev mafft

COPY . /make_prg
WORKDIR /make_prg
RUN python3 -m pip install --no-cache-dir .
# TEST
RUN python3 -m pip install --no-cache-dir nose hypothesis pytest
RUN PYTHONPATH="." pytest ./tests/ && make_prg --help
# CLEANUP
WORKDIR /
RUN python3 -m pip uninstall -y nose hypothesis pytest && \
    rm -rf /make_prg

# so that the use of python uses the python3 in the container
RUN apt install -y python3-venv
RUN python3 -m venv /venv
ENV PATH=/venv/bin:$PATH

RUN apt clean

CMD ["make_prg"]
