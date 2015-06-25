#simple bash function to assign command line arguments to variables
#example usage:
#  source /path/to/assign_args.sh
#  ARGS='|--option1|--option2|--option3|'
#  option3=DEFAULT_VALUE
#  assign_args ${ARGS} "$@"

assign_args()
{
    set -ue
    set -o pipefail
    
    CMD_LINE_OPTIONS=$1
    shift

    while [[ $# > 0 ]];do
        key=$1
        if echo ${CMD_LINE_OPTIONS} | grep -q -e "|${key}|"; then
            eval ${key:2}="$2"
        else
            echo "Command line options are: ${CMD_LINE_OPTIONS}" | tr '|' ' '
            exit 1
        fi
        shift;shift
    done
    
    TEST_ECHO=$(echo ${CMD_LINE_OPTIONS} | sed 's/|--/ \$/g' | tr '|' ' ')
    eval echo $TEST_ECHO > /dev/null
}

