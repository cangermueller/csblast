/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
