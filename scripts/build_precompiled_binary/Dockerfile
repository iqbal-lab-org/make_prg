# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/make_prg_precompiled_binary_builder:0.4.0
# This container has the python environment to build a static binary for make_prg
# To build: docker build . -t leandroishilima/make_prg_precompiled_binary_builder:0.4.0
FROM ubuntu:18.04
RUN apt update && apt install -y python3.7 python3.7-dev python3-pip libbz2-dev liblzma-dev graphviz graphviz-dev
COPY ./requirements.txt /requirements.txt
RUN python3.7 -m pip install --upgrade wheel pip && python3.7 -m pip install --prefer-binary -r requirements.txt
RUN python3.7 -m pip install pyinstaller
RUN cp -vr /usr/bin/* /usr/sbin/
