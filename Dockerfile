FROM alpine
MAINTAINER Michael Hall <michael@mbh.sh>

# DEPENDENCIES
RUN apk update && apk upgrade && \
    apk --no-cache add bash fd 'py3-scikit-learn>=0.19.1' 'py3-numpy>=1.14.0' py3-pip && \
    apk --no-cache add --virtual .builddeps gcc musl-dev
# INSTALL
COPY . /make_prg
WORKDIR /make_prg
RUN python3 -m pip install --no-cache-dir .
# TEST
RUN python3 -m pip install --no-cache-dir nose hypothesis pytest
RUN nosetests -v tests/ && make_prg --help
# CLEANUP
WORKDIR /
RUN python3 -m pip uninstall -y nose hypothesis pytest && \
    apk del .builddeps && \
    rm -rf /root/.cache && \
    rm -rf /make_prg

# so that the use of python uses the python3 in the container
RUN python3 -m venv /venv
ENV PATH=/venv/bin:$PATH

CMD ["python3", "-m", "make_prg"]