FROM ghcr.io/zhilongwei/wavein-petsc:3.25.4-r2-debug

ARG DEBIAN_FRONTEND=noninteractive
ARG USERNAME=ubuntu

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        openssh-client \
        catch2 \
        clang \
        clang-format \
        clang-tidy \
        clang-tools \
        clangd \
        gdb \
        lldb \
        sudo \
    && echo "${USERNAME} ALL=(root) NOPASSWD:ALL" > "/etc/sudoers.d/${USERNAME}" \
    && chmod 0440 "/etc/sudoers.d/${USERNAME}" \
    && rm -rf /var/lib/apt/lists/*

ENV CC=gcc \
    CXX=g++

USER ${USERNAME}
WORKDIR /workspaces/WaveIn
