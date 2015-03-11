if [ $# -eq 0 ]; then
    echo "Usage: $0 file.ddsg"
    exit
fi

if [ ! -f $1 ]; then
    echo "input file doesn't exist $1"
    exit
fi

if [ ! -f $1.hcn ]; then
    echo "Node ordering not found. Run first ./run-node-order.sh $1"
    exit
fi
./release/contraction -c -f $1 -h $1.hcn -C $1.ch
