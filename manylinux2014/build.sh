#!/bin/bash -ve

HERE=`dirname "$0"`
cd $HERE/..

for i in 36 37 38 39; do
    docker build -f manylinux2014/wheel${i}.docker . -t rmextract${i}
    docker run -v `pwd`/manylinux2014:/manylinux2014 rmextract${i} sh -c "cp /output/*.whl /manylinux2014/."
done

