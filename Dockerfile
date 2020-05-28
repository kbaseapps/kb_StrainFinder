FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update
RUN pip install openopt numpy scipy FuncDesigner DerApproximator Cython


# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module
RUN make all

RUN git clone https://bitbucket.org/yonatanf/strainfinder
WORKDIR /kb/module/strainfinder
RUN python setup_cython.py build_ext --inplace
ENV PATH $PATH:/kb/module/strainfinder

WORKDIR /kb/module

RUN git clone https://github.com/vcftools/vcftools.git
WORKDIR /kb/module/vcftools
RUN ./autogen.sh
RUN ./configure
RUN make
RUN make install

WORKDIR /kb/module
ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
