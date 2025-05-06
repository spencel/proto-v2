FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN pip install --upgrade pip

RUN pip install psycopg2-binary

WORKDIR /home/docker_postgresql
COPY psycopg.py psycopg.py 
COPY db_loader.py db_loader.py

CMD python3 db_loader.py