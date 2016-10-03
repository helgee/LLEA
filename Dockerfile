FROM centos:7

RUN yum update -y && yum install -y cmake gcc-gfortran
COPY . /llea
RUN mkdir /llea/build && mkdir /llea/debug
RUN cd /llea/build && cmake -DCMAKE_BUILD_TYPE=Release ..
RUN cd /llea/debug && cmake -DCMAKE_BUILD_TYPE=Debug ..
RUN make -C /llea/build
RUN make -C /llea/debug
RUN make -C /llea/build test
RUN make -C /llea/debug test
