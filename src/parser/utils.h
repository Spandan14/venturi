#pragma once

#include <fstream>
#include <sstream>

const std::string VENTURI_LANG_GRAMMAR_PATH = "../../../lang/venturi_lang.peg";

inline std::string load_text_file(const std::string &path) {
  std::ifstream ifs(path);
  if (!ifs)
    throw std::runtime_error("Cannot open: " + path);
  std::stringstream ss;
  ss << ifs.rdbuf();
  return ss.str();
}
