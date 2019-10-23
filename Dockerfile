FROM python:3.6.9-buster

RUN apt-get update && apt-get install -y python3 python3-pip

WORKDIR /usr/share/git/system_reliability
COPY . .

RUN pip3 install -r requirements.txt
