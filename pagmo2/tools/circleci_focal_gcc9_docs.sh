#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -y -q -p $deps_dir c-compiler cxx-compiler cmake eigen nlopt ipopt boost-cpp tbb tbb-devel python=3.8 sphinx'<4.1.0' sphinx_rtd_theme breathe doxygen graphviz
source activate $deps_dir

# Create the build dir and cd into it.
mkdir build
cd build

# GCC build.
cmake ../ -DCMAKE_BUILD_TYPE=Release -DPAGMO_BUILD_TESTS=yes -DPAGMO_WITH_EIGEN3=yes -DPAGMO_WITH_NLOPT=yes -DPAGMO_WITH_IPOPT=yes -DPAGMO_ENABLE_IPO=yes
make -j2 VERBOSE=1
ctest -V -j4

# Build the documentation.

# Doxygen.
cd ../doc/doxygen
export DOXYGEN_OUTPUT=`doxygen 2>&1 >/dev/null`;
if [[ "${DOXYGEN_OUTPUT}" != "" ]]; then
    echo "Doxygen encountered some problem:";
    echo "${DOXYGEN_OUTPUT}";
    exit 1;
fi
echo "Doxygen ran successfully";

# Copy the images into the xml output dir (this is needed by sphinx).
cp images/* xml/;
cd ../sphinx/;
make html linkcheck

if [[ ! -z "${CI_PULL_REQUEST}" ]]; then
    echo "Testing a pull request, the generated documentation will not be uploaded.";
    exit 0;
fi

if [[ "${CIRCLE_BRANCH}" != "master" ]]; then
    echo "Branch is not master, the generated documentation will not be uploaded.";
    exit 0;
fi

# Check out the gh_pages branch in a separate dir.
cd ../..
git config --global push.default simple
git config --global user.name "CircleCI"
git config --global user.email "bluescarni@gmail.com"
set +x
git clone "https://${GH_TOKEN}@github.com/esa/pagmo2.git" pagmo2_gh_pages -q
set -x
cd pagmo2_gh_pages
git checkout -b gh-pages --track origin/gh-pages;
git rm -fr *;
mv ../doc/sphinx/_build/html/* .;
git add *;
# We assume here that a failure in commit means that there's nothing
# to commit.
git commit -m "Update Sphinx documentation, commit ${CIRCLE_SHA1} [skip ci]." || exit 0
PUSH_COUNTER=0
until git push -q
do
    git pull -q
    PUSH_COUNTER=$((PUSH_COUNTER + 1))
    if [ "$PUSH_COUNTER" -gt 3 ]; then
        echo "Push failed, aborting.";
        exit 1;
    fi
done

set +e
set +x
