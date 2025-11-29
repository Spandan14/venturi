#pragma once

#include "peglib.h"
#include <engine/sim2d.h>
#include <fluxlang/utils.h>

class Parser {
public:
  Parser() = default;
  ~Parser() = default;

  void parse(const char *filename, Simulation2D &sim);
  // void parse(const std::string &src, Simulation3D &sim);

private:
  static peg::parser parser;

  static peg::parser _load_grammar(const std::string &grammar_path);
};
