FROM amazonlinux:latest

LABEL maintainer="nuno.fonseca at gmail.com"
RUN yum update -y && yum install -y git gcc wget make R && yum clean all
RUN git clone https://github.com/nunofonseca/sv2gf.git && cd sv2gf && ./run_tests.sh
