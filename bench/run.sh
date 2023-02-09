SCRIPT_DIR=$(dirname "$0")
RUN_BIN=$SCRIPT_DIR/bench.exe
clang -g -Wall -Wextra -Wno-missing-field-initializers $SCRIPT_DIR/bench.c -o $RUN_BIN -lpthread && $RUN_BIN $@