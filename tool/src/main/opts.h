#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C"
{
#include <getopt.h>
}
#include "util/exception.h"

struct options
{
    char input[255];
    char inputangle[255];
    char txtinput[255]; // txt输入文件
    float rotation_angle;
    float diameter;
    int verbose;
    bool testmode;
};

inline int GetOpts(int argc, char **argv, options *opts_)
{
    static struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"input", required_argument, NULL, 'i'},
        {"initangle", required_argument, NULL, 'a'},
        {"diameter", required_argument, NULL, 'd'},
        {"rotationangle", required_argument, NULL, 'r'},
        {"verbose", required_argument, NULL, 'v'},
        {"testmode", required_argument, NULL, 't'},
        {"txtinput", required_argument, NULL, 'x'},
        {NULL, 0, NULL, 0}};

    // 更新参数数量检查(删掉 -n 和 -o 后更少参数)
    if ((argc < 3) || (argc > 14))
    {
        EX_TRACE("[-i INPUT][-a INITANGLE][-d DIAMETER][-x TXTINPUT]([-r ROT])([-v VERBOSE])([-t TESTMODE])\n");
        return -1;
    }

    int ch;
    while ((ch = getopt_long(argc, argv, "hti:i:a:d:r:v:x:", longopts, NULL)) != -1)
    {
        switch (ch)
        {
        case '?':
            EX_TRACE("Invalid option '%s'.", argv[optind - 1]);
            return -1;

        case ':':
            EX_TRACE("Missing option argument for '%s'.", argv[optind - 1]);
            return -1;

        case 'h':
            EX_TRACE("  [--input(-i) INPUT]\n"
                     "  [--initangle(-a) INIT ANGLE]\n"
                     "  [--diameter(-d) MARKER DIAMETER]\n"
                     "  [--txtinput(-x) TXT INPUT]\n"
                     "  ([--rotationangle(-r) ROT])\n"
                     "  ([--verbose(-v) 0|1])\n"
                     "  ([--testmode(-t)])\n");
            return 0;

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if (iss.fail())
            {
                EX_TRACE("Invalid -i argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'a':
        {
            std::istringstream iss(optarg);
            iss >> opts_->inputangle;
            if (iss.fail())
            {
                EX_TRACE("Invalid -a argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'd':
        {
            std::istringstream iss(optarg);
            iss >> opts_->diameter;
            if (iss.fail())
            {
                EX_TRACE("Invalid -d argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'x':
        {
            std::istringstream iss(optarg);
            iss >> opts_->txtinput;
            if (iss.fail())
            {
                EX_TRACE("Invalid -x argument '%s'.", optarg);
                return -1;
            }
        }
        break;

        case 'r':
        {
            std::istringstream iss(optarg);
            iss >> opts_->rotation_angle;
            if (iss.fail())
                EX_TRACE("Invalid -r argument '%s'.", optarg);
        }
        break;

        case 'v':
        {
            std::istringstream iss(optarg);
            iss >> opts_->verbose;
            if (iss.fail())
                EX_TRACE("Invalid -v argument '%s'.", optarg);
        }
        break;

        case 't':
            opts_->testmode = true;
            break;

        default:
            assert(false);
        }
    }
    return 1;
}

#endif
