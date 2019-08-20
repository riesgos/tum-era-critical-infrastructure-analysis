FROM python:3.6.9-buster

RUN apt-get update && apt-get install -y python3 python3-pip && \
    pip3 install numpy networkx

WORKDIR /usr/share/git/system_reliability
COPY . .
