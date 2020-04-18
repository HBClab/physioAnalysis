#!/bin/sh

set -e

# Generate Dockerfile.
generate_docker() {
  docker run --rm kaczmarj/neurodocker:0.5.0 generate docker \
    --base=ubuntu:xenial-20161213 \
    --pkg-manager=apt \
    --install net-tools git \
    --user=neuro \
    --workdir="/home/neuro" \
    --run 'mkdir /home/neuro/physioAnalysis' \
    --copy . /home/neuro/physioAnalysis \
    --miniconda create_env='physio' \
                yaml_file='/home/neuro/physioAnalysis/environment.yml' \
    --env "SHELL=/bin/bash" \
    --run "curl -o /tmp/code-server.tar.gz -SL https://github.com/cdr/code-server/releases/download/3.0.2/code-server-3.0.2-linux-x86_64.tar.gz" \
    --run "mkdir -p /opt/codeserver && tar -xvf /tmp/code-server.tar.gz -C /opt/codeserver --strip-components=1" \
    --run '/opt/codeserver/code-server --install-extension eamodio.gitlens && /opt/codeserver/code-server --install-extension ms-python.python' \
    --expose 8080 \
    # --entrypoint '/opt/codeserver/code-server --auth none --host 0.0.0.0 /home/neuro/physioAnalysis'
}

generate_docker > Dockerfile

docker build -t jdkent/physio-analysis:dev .
# singularity build hbclab_nibetaseries.simg Singularity
