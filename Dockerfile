FROM containers.mathworks.com/matlab-runtime:r2023a

ENV DEBIAN_FRONTEND=noninteractive
ENV MCR_CACHE_ROOT=/tmp/mcr_cache
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility

ENV AGREE_TO_MATLAB_RUNTIME_LICENSE=yes

ENV MCR_ROOT=/opt/matlabruntime/R2023a
ENV LD_LIBRARY_PATH=${MCR_ROOT}/runtime/glnxa64:${MCR_ROOT}/bin/glnxa64:${MCR_ROOT}/sys/os/glnxa64:${MCR_ROOT}/extern/bin/glnxa64
ENV XAPPLRESDIR=${MCR_ROOT}/X11/app-defaults

# Install AWS CLI and dependencies
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    groff \
    less \
    && rm -rf /var/lib/apt/lists/*

# Install AWS CLI v2
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf aws awscliv2.zip

COPY bin/ /app/
COPY data/ /app/data/

RUN mkdir -p /app/data/solutions

RUN chmod +x /app/*

WORKDIR /app

ENTRYPOINT ["/app/run_simulator.sh", "/opt/matlabruntime/R2023a"]

CMD ["MARIE_runner", "data/inputs/inp_Duke_StadiumTriangular.json"]

# CMD ["MARIE_runner", "data/inputs/inp_Duke_StadiumTriangular.json", "s3-bucket-name"]