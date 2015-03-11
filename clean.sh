if [ $# -eq 0 ]; then
    echo "deletes all generated files"
    echo "Usage: $0 file.ddsg"
    exit
fi

rm -f $1.sgr $1.*.log $1.hcn -C $1.ch $1.test-lengths
