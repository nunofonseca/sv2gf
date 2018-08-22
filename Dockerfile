FROM ubuntu:16.04

LABEL maintainer="nuno.fonseca at gmail.com"
RUN apt-get update && apt-get install -y git gcc wget make r-base && rm -rf /var/lib/apt/lists/
#RUN yum update -y && yum install -y git gcc wget make R-base && yum clean all
RUN git clone https://github.com/nunofonseca/sv2gf.git && cd sv2gf && ./run_tests.sh
