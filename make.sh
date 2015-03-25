MODE=debug
if [ ${#@} -eq 1 ]; then
    MODE=$1
fi

if [ $MODE == "release" ]; then
    OPT_RELEASE_ARGS="-DCMAKE_BUILD_TYPE=Release"
elif [ $MODE == "debug" ]; then
    :
else
    echo "usage: ./make.sh [debug|release]"
    exit
fi

echo Compiling project in $MODE mode
rm -rf $MODE
mkdir $MODE
cd $MODE
cmake $OPT_RELEASE_ARGS ..
make
cd ..
