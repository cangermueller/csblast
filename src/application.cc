// Copyright 2009, Andreas Biegert

#include "application.h"

#include <cstdio>
#include <cstdlib>

#include <string>

#include "globals.h"
#include "exception.h"
#include "getopt_pp.h"
#include "log.h"

using std::string;

namespace cs {

Application* Application::instance_;

const char* Application::kVersionNumber = "2.0.1";

const char* Application::kCopyright =
  "Copyright (c) 2008 Andreas Biegert, Johannes Soding, and LMU Munich";

Application::Application()
    : log_level_(Log::to_string(Log::from_int(LOG_MAX_LEVEL))),
      log_file_fp_(NULL) {
  // Register the app. instance
  if (instance_)
    throw Exception("Second instance of Application is prohibited");
  instance_ = this;
}

Application::~Application() {
  if (log_file_fp_ && log_file_fp_ != stderr)
    fclose(log_file_fp_);
}

int Application::main(int argc,
                      char* argv[],
                      FILE* fout,
                      const string& name) {
  int status = 0;

  output_fp_     = fout;
  app_name_      = name;
  log_file_name_ = name + ".log";

  // Prepare command line parsing
  GetOpt_pp options(argc, argv, Include_Environment);
  options.exceptions_all();

  try {
    // Print usage?
    if (argc < 2 || argv[1][0] == '?' ||
        options >> OptionPresent(' ', "help")) {
      PrintHelp();
      return 1;
    }

    if (kDebug) {
      // Process logging options
      options >> Option(' ', "log-level", log_level_, log_level_);
      Log::reporting_level() = Log::from_string(log_level_);
      options >> Option(' ', "log-file", log_file_name_, log_file_name_);
      log_file_fp_ = fopen(log_file_name_.c_str(), "w");
      Log::stream() = log_file_fp_;
    }

    // Let subclasses parse the command line options
    ParseOptions(&options);

    // Run application
    status = Run();

  } catch(const std::exception& e) {
    LOG(ERROR) << e.what();
    fprintf(fout, "\n%s\n", e.what());
    return 1;
  }

  return status;
}

void Application::PrintHelp() const {
  fprintf(stream(), "%s version %s\n", app_name_.c_str(), kVersionNumber);
  PrintBanner();
  fprintf(stream(), "%s\n\n", kCopyright);
  PrintUsage();
  fputs("\nOptions:\n", stream());
  PrintOptions();

  if (kDebug) {
    fprintf(stream(), "  %-30s %s (def=%s)\n", "    --log-level <level>",
            "Maximal reporting level for logging", log_level_.c_str());
    fprintf(stream(), "  %-30s %s (def=%s)\n", "    --log-file <file>",
            "Output file for logging", log_file_name_.c_str());
  }
}

}  // namespace cs
