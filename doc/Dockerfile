FROM continuumio/anaconda3:latest

ENV SRC_HOME /usr/src/app

RUN mkdir -p $SRC_HOME
WORKDIR $SRC_HOME

COPY requirements.txt $SRC_HOME

RUN conda install --file=requirements.txt

ADD . $SRC_HOME
