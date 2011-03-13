// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "utils.h"

namespace cs {

// Reads all files in 'path' and pushes them onto given vector.
void GetAllFiles(const std::string& path,
                 std::vector<std::string>& files,
                 const std::string& ext) {
  DIR *dir = opendir(path.c_str());
  struct dirent *dp;  // returned from readdir()
  if (dir == NULL) return;  // error in opendir

  while ((dp = readdir (dir)) != NULL) {
    std::string filename(dp->d_name);
    // Leave out directories, '.', '..', and files with wrong extension
    if (filename[0] == '.' || IsDirectory(filename) ||
        (!ext.empty() && GetFileExt(filename) != ext)) continue;
    // Push filename onto vector
    files.push_back(std::string(dp->d_name));
  }
  closedir(dir);
  return;
}

}  // namespace cs
