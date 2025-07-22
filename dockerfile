FROM openjdk:21-slim

# Instala wget, unzip, curl, y libgomp1 (¡aquí está la clave!)
RUN apt-get update && apt-get install -y wget unzip git curl libgomp1

# (Resto del Dockerfile igual...)
COPY 2.7.11b.tar.gz /tmp/2.7.11b.tar.gz

RUN mkdir /opt/STAR && \
    tar -xzf /tmp/2.7.11b.tar.gz -C /opt/STAR && \
    ln -s /opt/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR /usr/local/bin/STAR

RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -O /usr/local/bin/picard.jar
RUN echo '#!/bin/bash\nexec java -jar /usr/local/bin/picard.jar "$@"' > /usr/local/bin/picard && chmod +x /usr/local/bin/picard

RUN wget https://github.com/broadinstitute/Drop-seq/releases/download/v3.0.2/dropseq-3.0.2.zip -O /tmp/dropseq-3.0.2.zip && \
    unzip /tmp/dropseq-3.0.2.zip -d /opt && \
    ln -s /opt/dropseq-3.0.2 /opt/dropseq

ENV PATH="/opt/dropseq:/usr/local/bin:${PATH}"

WORKDIR /data





