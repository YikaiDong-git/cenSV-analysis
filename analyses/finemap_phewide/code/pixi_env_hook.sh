: "${CENSV_WORKSPACE:?set CENSV_WORKSPACE (see README)}"
: "${HOME:?set HOME (see README)}"
export PATH="${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin"
export CONDA_SHLVL=1
export CONDA_PREFIX=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default
export PIXI_PROJECT_VERSION=0.1.0
export PIXI_IN_SHELL=1
export PIXI_PROJECT_MANIFEST=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/pixi.toml
export PIXI_EXE=${HOME}/.pixi/bin/pixi
export PIXI_PROJECT_ROOT=${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515
export PIXI_PROJECT_NAME=bayesian_finemap
export CONDA_DEFAULT_ENV=bayesian_finemap
export PIXI_ENVIRONMENT_NAME=default
export PIXI_ENVIRONMENT_PLATFORMS=linux-64
export PIXI_PROMPT='(bayesian_finemap) '
. ${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default/etc/conda/activate.d/activate-r-base.sh
. ${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default/etc/conda/activate.d/libglib_activate.sh
. ${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default/etc/conda/activate.d/libxml2_activate.sh
source ${CENSV_WORKSPACE}/_shared/tools/bayesian_finemap/v1_20260515/.pixi/envs/default/share/bash-completion/completions/*

pixi() {
    local first_arg="${1-}"

    "${PIXI_EXE-}" "$@" || return $?

    case "${first_arg-}" in
    add | a | remove | rm | install | i)
        eval "$("$PIXI_EXE" shell-hook --change-ps1 false)"
        hash -r
        ;;
    esac || :

    return 0
}

export PS1="(bayesian_finemap) ${PS1:-}"
