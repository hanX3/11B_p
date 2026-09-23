#ifndef OutputPath_H
#define OutputPath_H 1

#include <filesystem>
#include <string>

namespace HBOutputPath
{
  std::filesystem::path LocateProjectRoot();
  std::filesystem::path DataDirectory();
  std::filesystem::path MakeOutputFilePath(const std::string& prefix, const char* file_tag);
  void EnsureDataDirectory();
}

#endif
