#include <iostream>

#include "exception.h"
#include "log.h"
#include "util.h"

class Params
{
  public:
    Params() : num_states_(0)
    { }

    int num_states() { return num_states_; }
    Params& num_states(int n) { num_states_ = n; return *this; }

    void check()
    {
        if (num_states() == 0) throw cs::Exception("No value for number of HMM states provided!");
    }

  private:
    int num_states_;
};

std::ostream& usage(std::ostream &out);
void parse_options(const int argc, const char **argv, Params& params);
void check_params(Params& params);

int main(int argc, const char* argv[])
{
    if(argc < 2) {
        usage(std::cout);
        return 1;
    }

    try {
        Params params;
        parse_options(argc, argv, params);
        params.check();

        std::cout << "num_states=" << params.num_states() << std::endl;

    } catch(const std::exception& e) {
        LOG(ERROR) << e.what();
        std::cout << e.what() << std::endl;
        return 3;
    }
    return 0;
}

void parse_options(const int argc, const char **argv, Params& params)
{
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg == "-h" || arg == "-help" || arg == "--help") {
            usage(std::cout);
            exit(2);
        } else if (arg == "-K" || arg == "--num-states") {
            if (++i < argc) params.num_states(atoi(argv[i]));
            else throw cs::Exception("No value available for '%s'", arg.c_str());
        } else if (arg == "--log-max-level") {
            if (++i < argc) Log::reporting_level() = Log::from_string(argv[1]);
            else throw cs::Exception("No value available for '%s'", arg.c_str());
        }
    }
}

std::ostream& usage(std::ostream& out)
{
    out << "Train an HMM of context-states on a dataset of alignments.\n";
    out << "(C) Andreas Biegert, Johannes Soding, and Ludwig-Maximillians University Munich\n\n";

    out << "Usage: cstrain -i <alifile> -K <num_states> [options]\n\n";

    out << "Options:\n";
    out << cs::strprintf("  %-30s %s\n", "-K, --num-states <int>", "Number of states in the HMM to be trained");
    out << cs::strprintf("  %-30s %s (def=%i)\n", "    --log-max-level <int>", "Maximal reporting level for logging", Log::reporting_level());
}

