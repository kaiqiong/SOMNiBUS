THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export THISDIR
export PROJECTROOT=$THISDIR/../
export DOCKER_BUILDKIT=1

echo "THISDIR: " $THISDIR
echo "PROJECTROOT: " $PROJECTROOT
docker build $PROJECTROOT --file Dockerfile -t dk-somnibus
