if [ $# -eq 0 ]; then
    echo "Usage: $0 file.ddsg"
    exit
fi

if [ ! -f $1 ]; then
    echo "input file doesn't exist $1"
    exit
fi

./release/node-order -s -p -c -f $1 -o $1.hcn -l $1.node-order.log -Z $1.sgr -x 190 -y 1 -p 1000 -n 1000 -V 60 -S 145 -e 70
