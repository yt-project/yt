FROM fedora:latest
MAINTAINER Kacper Kowalik "xarthius.kk@gmail.com"

ENV DEBIAN_FRONTEND noninteractive

RUN yum update -y

# add a user
RUN useradd -D --shell=/bin/bash && \
    useradd -m ytuser && \
    echo "ytuser:secret" | chpasswd

WORKDIR /home/ytuser
RUN curl -O https://bitbucket.org/yt_analysis/yt/raw/tip/doc/install_script.sh && \
    sed -i -e '/^MAKE_PROCS=/ s/""/"-j4"/' install_script.sh

RUN yum install -y `grep -oP '(?<=yum install ).*?(?="|$) | tr "\n" " "' install_script.sh`

USER ytuser
RUN echo "" | bash install_script.sh
