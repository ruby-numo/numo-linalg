FROM andrewosh/binder-base

MAINTAINER Kozo Nishida <knishida@riken.jp>

USER root

RUN apt-get update
RUN apt-get install -y build-essential ruby ruby-dev rake git libzmq3 libzmq3-dev libgsl0-dev libtool autoconf automake zlib1g-dev libopenblas-dev && apt-get clean
RUN ln -s /usr/bin/libtoolize /usr/bin/libtool # See https://github.com/zeromq/libzmq/issues/1385
RUN git clone git://github.com/ruby-numo/narray; git clone git://github.com/ruby-numo/linalg; 
RUN gem update --no-document --system && gem install --no-document bundler daru iruby nyaplot pry rbczmq
RUN cd narray; gem build numo-narray.gemspec; gem install numo-narray-0.9.0.2.gem; cd ../linalg; rake build; gem install pkg/numo-linalg-0.0.1.gem -- --with-openblas

USER main

RUN iruby register
