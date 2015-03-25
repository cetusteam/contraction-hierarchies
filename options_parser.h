#pragma once

#include <string>
#include <iostream>
#include <getopt.h>

namespace {
using namespace std;
}

class options_parser_t {
public:

    void
    parse_options(int argc, char **argv) {
        int opt;

        while ((opt = getopt(argc, argv, opts_str)) != EOF) {
            handled = false;
            parse();

            if (!handled) {
                usage(argv[0]);
                exit(1);
            }
        }

        if (optind != argc) {
            printf("Invalid argument %s\n", argv[optind]);
            exit(1);
        }
        run_checks();
    }

protected:

    bool handled;
    const char *opts_str;
    options_parser_t(const char *opts_str) {
        this->opts_str = opts_str;
    }

    void
    set_option(char c, string *str) {
        if (optopt != c) {
            return;
        }
        *str = optarg; handled = true;
    }

    template<typename Fn>void
    if_option(char c, Fn f) {
        if (optopt != c) {
            return;
        }
        f(); handled = true;
    }

    void
    check_not_empty(const string& arg, const string& err_msg) {
        if (arg.empty()) {
            cerr << err_msg << endl;
            exit(1);
        }
    }

    virtual void parse()              = 0;
    virtual void run_checks()         = 0;
    virtual void usage(const char *p) = 0;
};
