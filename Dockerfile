FROM containers.mathworks.com/matlab-runtime:r2023a

ENV DEBIAN_FRONTEND=noninteractive
ENV MCR_CACHE_ROOT=/tmp/mcr_cache
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility

ENV AGREE_TO_MATLAB_RUNTIME_LICENSE=yes

ENV MCR_ROOT=/opt/matlabruntime/R2023a
ENV LD_LIBRARY_PATH=${MCR_ROOT}/runtime/glnxa64:${MCR_ROOT}/bin/glnxa64:${MCR_ROOT}/sys/os/glnxa64:${MCR_ROOT}/extern/bin/glnxa64
ENV XAPPLRESDIR=${MCR_ROOT}/X11/app-defaults

COPY bin/ /app/
COPY data/ /app/data/

RUN ln -sf /app/data/bodies /app/RHBM

RUN chmod +x /app/*

WORKDIR /app

ENTRYPOINT ["/app/run_simulator.sh", "/opt/matlabruntime/R2023a"]

CMD ["MARIE_runner", "data/inputs/inp_Duke_StadiumTriangular.json"]