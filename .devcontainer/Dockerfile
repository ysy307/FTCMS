# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# --- Environment variables and build arguments ---
ENV LANG=ja_JP.UTF-8 \
    LANGUAGE=ja_JP:ja \
    LC_ALL=ja_JP.UTF-8 \
    TZ=Asia/Tokyo

ENV SHELL=/usr/bin/zsh

ARG USERNAME=devuser
ARG USER_UID=1000
ARG USER_GID=${USER_UID}

# --- 1. Install base packages ---
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    wget gnupg htop lsb-release curl build-essential ca-certificates sudo vim git \
    gfortran python3 python3-pip python3-dev zsh sed software-properties-common \
    ninja-build graphviz linux-tools-common linux-tools-generic \
    libx11-dev libxext-dev libxrender-dev libxt-dev pkg-config \
    adwaita-icon-theme-full gnuplot valgrind kcachegrind gdb tree \
    openmpi-bin libopenmpi-dev \
    locales tzdata tmux flex libgsl-dev && \
    locale-gen ja_JP.UTF-8 && \
    update-locale LANG=ja_JP.UTF-8 LC_ALL=ja_JP.UTF-8 && \
    ln -snf /usr/share/zoneinfo/Asia/Tokyo /etc/localtime && \
    echo Asia/Tokyo > /etc/timezone && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# --- 2. Add external repositories ---
RUN wget -qO- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor > /usr/share/keyrings/oneapi-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" > /etc/apt/sources.list.d/oneAPI.list && \
    curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | gpg --dearmor > /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /" > /etc/apt/sources.list.d/nvhpc.list && \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor > /usr/share/keyrings/kitware-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" > /etc/apt/sources.list.d/kitware.list && \
    mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor > /etc/apt/keyrings/nodesource.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_22.x nodistro main" > /etc/apt/sources.list.d/nodesource.list

# --- 3. Install compilers, CMake, and Node.js ---
RUN apt-get update -y && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    intel-oneapi-compiler-fortran \
    intel-oneapi-hpc-toolkit \
    nvhpc-25-3 \
    cmake \
    nodejs && \
    node --version && \
    npm --version && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# --- 4. Environment setup scripts ---
RUN echo 'export NVARCH=$(uname -s)_$(uname -m)' > /etc/profile.d/nvhpc.sh && \
    echo 'export NVCOMPILERS=/opt/nvidia/hpc_sdk' >> /etc/profile.d/nvhpc.sh && \
    echo 'export PATH=$NVCOMPILERS/$NVARCH/25.3/compilers/bin:$PATH' >> /etc/profile.d/nvhpc.sh && \
    echo 'export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.3/compilers/man' >> /etc/profile.d/nvhpc.sh && \
    echo 'export PATH=$NVCOMPILERS/$NVARCH/25.3/comm_libs/mpi/bin:$PATH' >> /etc/profile.d/nvhpc.sh && \
    echo 'export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/25.3/comm_libs/mpi/man' >> /etc/profile.d/nvhpc.sh && \
    echo '[ -f /opt/intel/oneapi/setvars.sh ] && source /opt/intel/oneapi/setvars.sh' > /etc/profile.d/oneapi.sh

# --- 5. Install Python tools ---
RUN pip3 install --upgrade pip && \
    pip3 install --force-reinstall git+https://github.com/fortran-lang/fortls && \
    pip3 install --upgrade fprettify fypp ford

# --- 6. Create regular user ---
RUN groupadd --gid ${USER_GID} ${USERNAME} && \
    useradd --uid ${USER_UID} --gid ${USER_GID} -m -s ${SHELL} ${USERNAME} && \
    echo ${USERNAME} ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/${USERNAME} && \
    chmod 0440 /etc/sudoers.d/${USERNAME}

# --- 7. Configure Zsh and Prezto as devuser ---
USER ${USERNAME}
WORKDIR /home/${USERNAME}

RUN git clone --recursive https://github.com/sorin-ionescu/prezto.git "${ZDOTDIR:-$HOME}/.zprezto" && \
    git clone --recurse https://github.com/belak/prezto-contrib "${ZDOTDIR:-$HOME}/.zprezto/contrib"

RUN zsh -c 'setopt EXTENDED_GLOB && \
    for rcfile in "${ZDOTDIR:-$HOME}"/.zprezto/runcoms/^README.md(.N); do \
    ln -sf "$rcfile" "${ZDOTDIR:-$HOME}/.${rcfile:t}"; \
    done'

RUN sed -i "/'prompt'/c \\\ 'contrib-prompt' \\\\\n  'prompt'" .zpreztorc && \
    sed -i "s/theme 'sorin'/theme 'bart'/g" .zpreztorc

RUN echo 'export PATH="$HOME/.local/bin:$PATH"' >> .zshrc && \
    echo 'export PS1="%~ %# "' >> .zshrc && \
    echo "export FC=ifx" >> .zshrc && \
    echo "export F90=ifx" >> .zshrc && \
    echo "export F77=ifx" >> .zshrc && \
    echo "export CC=icx" >> .zshrc && \
    echo "export CXX=icpx" >> .zshrc && \
    echo "export MPIFC=mpiifx" >> .zshrc && \
    echo "export MPIF90=mpiifx" >> .zshrc && \
    echo "export MPIF77=mpiifx" >> .zshrc && \
    echo "export MPICC=mpiicc" >> .zshrc && \
    echo "export MPICXX=mpiicpc" >> .zshrc && \
    echo "export PETSC_DIR=/workspaces/FTCMS/third_party/.local" >> .zshrc && \
    echo "ulimit -s unlimited" >> .zshrc && \
    echo '[ -f /etc/profile.d/nvhpc.sh ] && source /etc/profile.d/nvhpc.sh' >> .zprofile && \
    echo '[ -f /etc/profile.d/oneapi.sh ] && source /etc/profile.d/oneapi.sh' >> .zprofile && \
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> .zprofile

# --- 8. Install AI CLIs ---
RUN mkdir -p "$HOME/.local" && \
    npm config set prefix "$HOME/.local" && \
    npm install -g \
    @anthropic-ai/claude-code@latest \
    @openai/codex@latest \
    @google/gemini-cli@latest && \
    curl -fsSL https://antigravity.google/cli/install.sh | bash && \
    "$HOME/.local/bin/claude" --version && \
    "$HOME/.local/bin/codex" --version && \
    "$HOME/.local/bin/gemini" --version && \
    test -x "$HOME/.local/bin/agy"

# --- 9. Add a helper script for project-local Codex MCP registration ---

USER root

RUN cat > /usr/local/bin/register-codex-mcp <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

cd "${1:-/workspaces/FTCMS}"

if ! command -v claude >/dev/null 2>&1; then
    echo "claude command was not found."
    exit 1
fi

if ! command -v codex >/dev/null 2>&1; then
    echo "codex command was not found."
    exit 1
fi

claude mcp add --scope project codex -- codex mcp-server
claude mcp list
EOF

RUN chmod 0755 /usr/local/bin/register-codex-mcp && \
    ln -sf /home/${USERNAME}/.local/bin/claude /usr/local/bin/claude && \
    ln -sf /home/${USERNAME}/.local/bin/codex /usr/local/bin/codex && \
    chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}

USER ${USERNAME}

# --- 10. Default command ---

CMD ["/usr/bin/zsh"]