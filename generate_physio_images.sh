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
    --user=root \
    --run 'curl -o /tmp/code-server.tar.gz -SL https://github.com/cdr/code-server/releases/download/2.1472-vsc1.38.1/code-server2.1472-vsc1.38.1-linux-x86_64.tar.gz' \
    --run 'mkdir -p /opt/codeserver && tar -xvf /tmp/code-server.tar.gz -C /opt/codeserver --strip-components=1' \
    --user=neuro \
    --run '/opt/codeserver/code-server --install-extension eamodio.gitlens && /opt/codeserver/code-server --install-extension ms-python.python' \
    --entrypoint '/opt/codeserver/code-server /home/neuro/physioAnalysis'

}

generate_docker > Dockerfile

docker build -t jdkent/physio-analysis:dev .
# singularity build hbclab_nibetaseries.simg Singularity
