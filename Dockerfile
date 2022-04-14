FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

WORKDIR /kb/module

RUN apt-get update
#RUN apt-get -y install build-essential libssl-dev
#RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py | python  # fails
#RUN pip install openopt numpy scipy FuncDesigner DerApproximator Cython  # fails


# NOTE: can't use pip anymore for packages needed by StrainFinder: openopt FuncDesigner DerApproximator (because python 2.7?)

RUN pip install numpy scipy Cython

#
# prep dependencies needed by openopt FuncDesigner and DerApproximator packages
#

# manually install setproctitle (for sortedcontainers)
WORKDIR /kb/module
ENV PACKAGE='setproctitle-1.2.2'
ENV REMOTE_URL='https://files.pythonhosted.org/packages/a1/7f/a1d4f4c7b66f0fc02f35dc5c85f45a8b4e4a7988357a29e61c14e725ef86/setproctitle-1.2.2.tar.gz'
RUN curl -o ${PACKAGE}.tar.gz ${REMOTE_URL} && \
  tar xfz ${PACKAGE}.tar.gz && \
  rm ${PACKAGE}.tar.gz && \
  cd ${PACKAGE} && \
  python setup.py install

# manually install sortedcontainers (for openopt)
WORKDIR /kb/module
ENV PACKAGE='sortedcontainers-2.4.0'
ENV REMOTE_URL='https://files.pythonhosted.org/packages/e8/c4/ba2f8066cceb6f23394729afe52f3bf7adec04bf9ed2c820b39e19299111/sortedcontainers-2.4.0.tar.gz'
RUN curl -o ${PACKAGE}.tar.gz ${REMOTE_URL} && \
  tar xfz ${PACKAGE}.tar.gz && \
  rm ${PACKAGE}.tar.gz && \
  cd ${PACKAGE} && \
  python setup.py install

#
# Below are StrainFinder_v1 packages
#

# manually install openopt
WORKDIR /kb/module
ENV PACKAGE='openopt-0.5629'
ENV REMOTE_URL='https://files.pythonhosted.org/packages/50/aa/214fbf85ecf1f5949114f5c2a15d0f96cf18d668259624bd80e24174e982/openopt-0.5629.tar.gz'
RUN curl -o ${PACKAGE}.tar.gz ${REMOTE_URL} && \
  tar xfz ${PACKAGE}.tar.gz && \
  rm ${PACKAGE}.tar.gz && \
  cd ${PACKAGE} && \
  python setup.py install

# manually install FuncDesigner
WORKDIR /kb/module
ENV PACKAGE='FuncDesigner-0.5629'
ENV REMOTE_URL='https://files.pythonhosted.org/packages/fa/32/c8c94bdf9da20cdedf521a34f7e5e9958a4f6420539bbe44bd5ff0e6c6fd/FuncDesigner-0.5629.tar.gz'
RUN curl -o ${PACKAGE}.tar.gz ${REMOTE_URL} && \
  tar xfz ${PACKAGE}.tar.gz && \
  rm ${PACKAGE}.tar.gz && \
  cd ${PACKAGE} && \
  python setup.py install

# manually install DerApproximator
WORKDIR /kb/module
ENV PACKAGE='DerApproximator-0.5201'
ENV REMOTE_URL='https://files.pythonhosted.org/packages/83/a5/2d6d62b1796e4617849c5c528b532dae7e25b4846872ad88fd91844135ea/DerApproximator-0.5201.tar.gz'
RUN curl -o ${PACKAGE}.tar.gz ${REMOTE_URL} && \
  tar xfz ${PACKAGE}.tar.gz && \
  rm ${PACKAGE}.tar.gz && \
  cd ${PACKAGE} && \
  python setup.py install


# return urllib3 to 1.22 which has ProtocolError, required by baseclient.py
RUN rm -rf /usr/lib/python2.7/dist-packages/urllib3
RUN pip install urllib3==1.22
ENV PYTHONPATH=${PYTHONPATH}:/usr/local/lib/python2.7/dist-packages

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module
RUN make all

#RUN git clone https://bitbucket.org/yonatanf/strainfinder
RUN curl -o yonatanf-strainfinder-2868046ef82b.zip https://bitbucket.org/yonatanf/strainfinder/get/2868046ef82b.zip && \
    unzip yonatanf-strainfinder-2868046ef82b.zip && \
    ln -s yonatanf-strainfinder-2868046ef82b strainfinder
WORKDIR /kb/module/strainfinder
RUN python setup_cython.py build_ext --inplace
ENV PATH $PATH:/kb/module/strainfinder

WORKDIR /kb/module

# vcftools run by meta_decoder subcall
#RUN git clone https://github.com/vcftools/vcftools.git
#WORKDIR /kb/module/vcftools
#RUN ./autogen.sh
#RUN ./configure
#RUN make
#RUN make install

WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
