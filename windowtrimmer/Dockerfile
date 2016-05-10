FROM python:2.7
MAINTAINER alejandro.barrera@duke.edu

ADD src /src
WORKDIR /src

RUN chmod +x /src/windowTrimmer.py
ENV PATH /src/:$PATH

CMD ["windowTrimmer.py"]