#include "OutputPath.hh"

#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

namespace
{
  bool HasProjectMarkers(const std::filesystem::path& path)
  {
    return std::filesystem::exists(path / "CMakeLists.txt") && std::filesystem::exists(path / "HB.cc");
  }

  std::filesystem::path ParseCMakeHomeDirectory(const std::filesystem::path& cache_path)
  {
    std::ifstream input(cache_path);
    std::string line;

    while(std::getline(input, line)) {
      const std::string key = "CMAKE_HOME_DIRECTORY:";
      if(line.rfind(key, 0) != 0) {
        continue;
      }

      const auto equals_pos = line.find('=');
      if(equals_pos == std::string::npos || equals_pos + 1 >= line.size()) {
        continue;
      }

      return std::filesystem::path(line.substr(equals_pos + 1));
    }

    return {};
  }

  std::filesystem::path AbsolutePathFromEnv(const char* name)
  {
    const char* value = std::getenv(name);
    if(value == nullptr || value[0] == '\0') {
      return {};
    }

    return std::filesystem::absolute(std::filesystem::path(value));
  }
}

namespace HBOutputPath
{
  std::filesystem::path LocateProjectRoot()
  {
    const auto env_project_root = AbsolutePathFromEnv("HB_PROJECT_ROOT");
    if(!env_project_root.empty()) {
      return env_project_root;
    }

    auto current = std::filesystem::absolute(std::filesystem::current_path());

    for(auto path = current; !path.empty(); path = path.parent_path()) {
      if(HasProjectMarkers(path)) {
        return path;
      }

      const auto cache_path = path / "CMakeCache.txt";
      if(std::filesystem::exists(cache_path)) {
        const auto cmake_home = ParseCMakeHomeDirectory(cache_path);
        if(!cmake_home.empty()) {
          return std::filesystem::absolute(cmake_home);
        }
      }

      if(path == path.root_path()) {
        break;
      }
    }

    return current;
  }

  std::filesystem::path DataDirectory()
  {
    const auto env_data_dir = AbsolutePathFromEnv("HB_DATA_DIR");
    if(!env_data_dir.empty()) {
      return env_data_dir;
    }

    return LocateProjectRoot() / "data";
  }

  std::filesystem::path MakeOutputFilePath(const char* file_tag)
  {
    EnsureDataDirectory();
    return DataDirectory() / (std::string(file_tag) + ".root");
  }

  void EnsureDataDirectory()
  {
    const auto data_directory = DataDirectory();
    std::error_code error;
    std::filesystem::create_directories(data_directory, error);
    if(error) {
      throw std::runtime_error("failed to create data directory: " + data_directory.string() + ": " + error.message());
    }
  }
}
