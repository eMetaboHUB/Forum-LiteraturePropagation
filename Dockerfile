FROM continuumio/miniconda3:latest

# SETTINGS
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get install -y libcairo2-dev libjpeg-dev libgif-dev build-essential

