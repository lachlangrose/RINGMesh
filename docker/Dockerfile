FROM ubuntu
WORKDIR /root/

ENV RINGMESH_PATH /root/RINGMesh
WORKDIR $RINGMESH_PATH

RUN apt-get update && apt-get install -y \
gcc \
g++ \
git \
cmake \
&& apt-get clean

RUN git clone https://github.com/lachlangrose/RINGMesh

RUN cd RINGMesh && \
mkdir build  && \
cd build && \

cmake .. && \

cd Release && \
make
