// Copyright 2011, Andreas Biegert

#include "cs.h"
#include "cslast.h"

using std::string;

namespace cs {

#ifdef _WIN32
const char* CSLast::kLastExec = "lastal.exe";
#else
const char* CSLast::kLastExec = "lastal";
#endif

const char* CSLast::kIgnoreOptions = "";


CSLast::CSLast(const string& dbfile,
               const LastPssm* pssm,
               const CSLastOptions& opts)
        : dbfile_(dbfile), pssm_(pssm), opts_(opts), exec_path_() {}

int CSLast::Run(FILE* fout) {
    int status = 0;
    string basename, pssmfile, resultsfile;

    try {
        // Create unique basename for sequence and checkpoint file
        char name_template[] = "/tmp/csblast_XXXXXX";
        const int captured_fd = mkstemp(name_template);
        if (!captured_fd) throw Exception("Unable to create unique filename!");
        basename    = name_template;
        pssmfile    = basename + ".pssm";
        resultsfile = basename + ".out";

        WritePssm(pssmfile);
        FILE* fres = fopen(resultsfile.c_str(),"w+");
        if (!fres) throw Exception("Unable to open file '%s'!", resultsfile.c_str());

        // Run LAST with provided options
        string command(ComposeCommandString(pssmfile));
        FILE* last_out = popen(command.c_str(), "r");
        if (!last_out) throw Exception("Error executing '%s'", command.c_str());

        // Read LAST output and print to stdout
        int c;
        while (!feof(last_out)) {
            if ((c = fgetc(last_out)) != EOF && fout) {
                fputc(c, fout);
                fflush(fout);
                fputc(c, fres);
                fflush(fres);
            }
        }
        status = pclose(last_out);

        // Cleanup temporary files
        if (!basename.empty()) remove(basename.c_str());
        if (!pssmfile.empty()) remove(pssmfile.c_str());
        if (!resultsfile.empty()) remove(resultsfile.c_str());

    } catch(const std::exception&) { // something went wrong => cleanup before rethrow
        // Cleanup temporary files
        if (!basename.empty()) remove(basename.c_str());
        if (!pssmfile.empty()) remove(pssmfile.c_str());
        if (!resultsfile.empty()) remove(resultsfile.c_str());
        throw;
    }
    return status;
}

void CSLast::WritePssm(string filepath) const {
    FILE* fout = fopen(filepath.c_str(), "w");
    if (!fout) throw Exception("Unable to write to file '%s'!", filepath.c_str());
    pssm_->Write(fout);
    // pssm_->Write(stdout);
    fclose(fout);
}

string CSLast::ComposeCommandString(string pssmfile) const {
    string rv(exec_path_ + kLastExec);
    if (!IsRegularFile(rv))
        throw Exception("No lastal binary in directory '%s'!", exec_path_.c_str());

    rv = rv + " -Q5 ";

    typedef CSLastOptions::const_iterator OptsIter;
    for (OptsIter it = opts_.begin(); it != opts_.end(); ++it) {
        if (!strchr(kIgnoreOptions, it->first))
            rv = rv + " -" + it->first + " " + it->second;
    }
    rv = rv + " '" + dbfile_ + "' " + pssmfile;

    return rv;
}

}  // namespace cs
