THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export THISDIR
export PROJECTROOT=$THISDIR/../
export DOCKER_BUILDKIT=1

echo "THISDIR: " $THISDIR
echo "PROJECTROOT: " $PROJECTROOT
docker build $PROJECTROOT --file Dockerfile -t dk-somnibus:v1
docker tag dk-somnibus happyregistry.azurecr.io/methods/wes/dk-somnibus:v1
docker push happyregistry.azurecr.io/methods/wes/dk-somnibus:v1
