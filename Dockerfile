FROM quay.io/pypa/manylinux2014_x86_64
COPY scripts/fossil-linux-x64-2.20.tar.gz /fossil-linux-x64-2.20.tar.gz
RUN tar -xf /fossil-linux-x64-2.20.tar.gz; mv fossil /bin; rm fossil-linux-x64-2.20.tar.gz
WORKDIR /project
RUN python3.9 -m pip install poetry
RUN python3.10 -m pip install poetry
RUN python3.11 -m pip install poetry
RUN python3.12 -m pip install poetry
COPY pipeline.fossil /
RUN fossil open /pipeline.fossil 
CMD ["bash", "scripts/build_python3.10.sh"]
