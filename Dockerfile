FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update
RUN pip install openopt numpy scipy FuncDesigner DerApproximator Cython
#RUN apt-get -y install build-essential libssl-dev
 

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
