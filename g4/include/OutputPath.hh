#ifndef OutputPath_H
#define OutputPath_H 1

#include <filesystem>

namespace HBOutputPath
{
  std::filesystem::path LocateProjectRoot();
  std::filesystem::path DataDirectory();
  std::filesystem::path MakeOutputFilePath(const char* file_tag);
  void EnsureDataDirectory();
}

#endif
