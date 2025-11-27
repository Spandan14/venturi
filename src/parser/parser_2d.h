#pragma once

#include "peglib.h"
#include <engine/sim2d.h>
#include <parser/utils.h>

class Parser2D {
public:
  Parser2D();
  ~Parser2D() = default;

  void parse(const char *filename, Simulation2D &sim);

private:
  peg::parser _load_grammar(const std::string &grammar_path);
  peg::parser parser;
};
