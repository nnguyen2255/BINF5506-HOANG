FROM ubuntu:latest

RUN apt update

RUN apt upgrade -y

RUN apt install tree


CMD tree